// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  udm.c	                                                                                                    */
/*																			*/
/* Purpose:   user data manager                                                 */
/*																			*/
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   December 9, 1996                                                                          */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>

#include "general.h"
#include "debug.h"
#include "gm.h"
#include "ugenv.h"
#include "devices.h"
#include "np.h"

#include "udm.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAX_NAMES 99

#define READ_DR_VEC_FLAG(p,vt,i)        READ_FLAG((p)->data_status.VecReserv[vt][(i)/32],1<<((i)%32))
#define READ_DR_MAT_FLAG(p,vt,i)        READ_FLAG((p)->data_status.MatReserv[vt][(i)/32],1<<((i)%32))
#define SET_DR_VEC_FLAG(p,vt,i)         SET_FLAG((p)->data_status.VecReserv[vt][(i)/32],1<<((i)%32))
#define SET_DR_MAT_FLAG(p,vt,i)         SET_FLAG((p)->data_status.MatReserv[vt][(i)/32],1<<((i)%32))
#define CLEAR_DR_VEC_FLAG(p,vt,i)       CLEAR_FLAG((p)->data_status.VecReserv[vt][(i)/32],1<<((i)%32))
#define CLEAR_DR_MAT_FLAG(p,vt,i)       CLEAR_FLAG((p)->data_status.MatReserv[vt][(i)/32],1<<((i)%32))

#define A_REASONABLE_NUMBER                     100

/* vm decriptor lock status */
#define VM_LOCKED(p)               ((p)->locked)
#define VM_IS_UNLOCKED                          0
#define VM_IS_LOCKED                            1

/* for SwapPartInterfaceData */
#define SWAP_VEC_DATA(v,pf,pt)          {tmp = VVALUE(v,*pf); VVALUE(v,*pf) = VVALUE(v,*pt); VVALUE(v,*pt) = tmp;}
#define SWAP_MAT_DATA(m,pf,pt)          {tmp = MVALUE(m,*pf); MVALUE(m,*pf) = MVALUE(m,*pt); MVALUE(m,*pt) = tmp;}

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static char NoVecNames[MAX_VEC_COMP];
static char NoMatNames[2*MAX_MAT_COMP];

static INT VectorDirID;
static INT MatrixDirID;
static INT VectorVarID;
static INT MatrixVarID;

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   GetUniqueOTypeOfVType - get uniqe object type for vtype

   SYNOPSIS:
   INT GetUniqueOTypeOfVType (const FORMAT *fmt, INT vtype)

   PARAMETERS:
   .  fmt - FORMAT
   .  vtype - check this vtype

   DESCRIPTION:
   This function checks whether vtype uses exactly one object type
   and if so returns its ID.

   RETURN VALUE:
   INT
   .n    number of object if unique
   .n    -1 else
   D*/
/****************************************************************************/

INT GetUniqueOTypeOfVType (const FORMAT *fmt, INT vtype)
{
  INT i,found,objs,obj;

  objs = FMT_T2O(fmt,vtype);
  found = 0;
  for (i=0; i<MAXVOBJECTS; i++)
    if (objs & (1<<i))
    {
      found++;
      obj = i;
    }

  if (found==1)
    return (obj);
  else
    REP_ERR_RETURN (-1);
}

/****************************************************************************/
/*D
   GetUniquePartOfVType - get uniqe part for vtype

   SYNOPSIS:
   INT GetUniquePartOfVType (const MULTIGRID *mg, INT vtype)

   PARAMETERS:
   .  mg - multigrid
   .  vtype - check this vtype

   DESCRIPTION:
   This function checks whether vtype uses exactly one part
   and if so returns its ID.

   RETURN VALUE:
   INT
   .n    number of object if unique
   .n    -1 else
   D*/
/****************************************************************************/

INT GetUniquePartOfVType (const MULTIGRID *mg, INT vtype)
{
  FORMAT *fmt=MGFORMAT(mg);
  INT i,n,found,parts,part;

  n = BVPD_NPARTS(MG_BVPD(mg));
  parts = FMT_T2P(fmt,vtype);
  found = 0;
  for (i=0; i<n; i++)
    if (parts & (1<<i))
    {
      found++;
      part = i;
    }

  if (found==1)
    return (part);
  else
    REP_ERR_RETURN (-1);
}

/****************************************************************************/
/*D
   IsVDdefinedInAllObjects - check whether descriptor covers objects in all parts

   SYNOPSIS:
   INT IsVDdefinedInAllObjects (const MULTIGRID *mg, const VECDATA_DESC *vd, INT obj_flags)

   PARAMETERS:
   .  mg		 - multigrid
   .  vd		 - check this descriptor ...
   .  obj_flags - ... for these objects (bitwise coded)

   DESCRIPTION:
   This function checks whether a vector descriptor covers certain object types
   in all parts (i.e. in the whole domain)

   RETURN VALUE:
   INT
   .n    YES
   .n    NO
   D*/
/****************************************************************************/

INT IsVDdefinedInAllObjects (const MULTIGRID *mg, const VECDATA_DESC *vd, INT obj_flags)
{
  FORMAT *fmt=MGFORMAT(mg);
  INT tp,i,n,parts;

  parts = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
      if (FMT_T2O(fmt,tp) & obj_flags)
        parts |= FMT_T2P(fmt,tp);

  n = BVPD_NPARTS(MG_BVPD(mg));
  for (i=0; i<n; i++)
    if (!(parts & (1<<i)))
      /* not all parts covered */
      return (NO);

  return (YES);
}

/****************************************************************************/
/*D
   FillCompsForOType - fill a vtype component vector for a certain object type

   SYNOPSIS:
   INT FillCompsForOType (const FORMAT *fmt, INT otype, INT n, SHORT cmps[])

   PARAMETERS:
   .  otype - want dofs in this object type
   .  n - number of dofs
   .  cmps - resulting component vector

   DESCRIPTION:
   This function fills a vtype component vector with n components in each vtype
   using objects of 'otype'. The 'cmps' vector can be used for dynamic allocation
   of a VECDATA_DESC.

   EXAMPLE:
   .vb
   SHORT NComps[NVECTYPES];
   VECDATA_DESC *vd;

   FillCompsForOType(MGFORMAT(theMG),NODEVEC,1,NComps);
   if (AllocVDfromNCmp(theMG,fl,tl,NComps,NULL,&vd)) return(error);
   .ve

   RETURN VALUE:
   INT
   .n    number of object if unique
   .n    -1 else
   D*/
/****************************************************************************/

INT FillCompsForOType (const FORMAT *fmt, INT otype, INT n, SHORT cmps[])
{
  INT tp,otp;

  otp = 1<<otype;
  for (tp=0; tp<NVECTYPES; tp++)
    cmps[tp] = (FMT_T2O(fmt,tp) & otp) ? n : 0;

  return (0);
}

/****************************************************************************/
/****************************************************************************/
/*			here follows vector stuff										*/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*D
   ConstructVecOffsets - Calculate offsets for VEC_SCALARs

   SYNOPSIS:
   INT ConstructVecOffsets (const SHORT *NCmpInType, SHORT *offset);

   PARAMETERS:
   .  NCmpInType - the numbers of components is used for the calculation of offsets
   .  offset - array of length NVECOFFSETS (!)

   DESCRIPTION:
   This function calculates offsets in DOUBLE vector called 'VEC_SCALAR'.
   It describes the number of components of each abstract type.

   .n      offset[0] = 0
   .n      offset[1] - offset[0] number of components in first type
   .n      offset[2] - offset[1] number of components in second type
   .n      etc.

   All components of a vector data descriptor are mapped uniquely to
   one of the DOUBLE values in the 'VEC_SCALAR'.

   'VD_CMP_OF_TYPE(vd,type,n)'

   RETURN VALUE:
   INT
   .n    NUM_OK
   D*/
/****************************************************************************/

INT ConstructVecOffsets (const SHORT *NCmpInType, SHORT *offset)
{
  INT type;

  offset[0] = 0;
  for (type=0; type<NVECTYPES; type++)
    offset[type+1] = offset[type] + NCmpInType[type];

  return (NUM_OK);
}

/****************************************************************************/
/*D
   SetScalVecSettings - fill the scalar settings components of a VECDATA_DESC

   SYNOPSIS:
   INT SetScalVecSettings (VECDATA_DESC *vd)

   PARAMETERS:
   .  vd - fill this scalar settings components

   DESCRIPTION:
   This function fills the scalar settings components of a VECDATA_DESC.

   RETURN VALUE:
   INT
   .n    NUM_OK
   D*/
/****************************************************************************/

static INT SetScalVecSettings (VECDATA_DESC *vd)
{
  INT tp;

  VD_IS_SCALAR(vd) = FALSE;

  /* check number of components per type */
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
    {
      if (VD_NCMPS_IN_TYPE(vd,tp)!=1)
        return (NUM_OK);                                                        /* no scalar */
      else
        VD_SCALCMP(vd) = VD_CMP_OF_TYPE(vd,tp,0);
    }

  /* check location of components per type */
  VD_SCALTYPEMASK(vd) = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
    {
      VD_SCALTYPEMASK(vd) |= 1<<tp;
      if (VD_SCALCMP(vd)!=VD_CMP_OF_TYPE(vd,tp,0))
        return (NUM_OK);                                                        /* no scalar */
    }

  VD_IS_SCALAR(vd) = TRUE;
  return (NUM_OK);
}

static INT SetCompactTypesOfVec (VECDATA_DESC *vd)
{
  FORMAT *fmt;
  INT tp;

  /* fill bitwise fields */
  fmt = MGFORMAT(VD_MG(vd));
  VD_DATA_TYPES(vd) = VD_OBJ_USED(vd) = 0;
  VD_MAXTYPE(vd) = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
    {
      VD_MAXTYPE(vd) = tp;
      VD_DATA_TYPES(vd) |= BITWISE_TYPE(tp);
      VD_OBJ_USED(vd)   |= FMT_T2O(fmt,tp);
    }
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
      break;
  VD_MINTYPE(vd) = tp;

  return (0);
}

static INT VDCompsSubsequent (const VECDATA_DESC *vd)
{
  INT tp,i;

  for (tp=0; tp<NVECTYPES; tp++)
    for (i=0; i<VD_NCMPS_IN_TYPE(vd,tp); i++)
      if (VD_CMP_OF_TYPE(vd,tp,i)!=VD_CMP_OF_TYPE(vd,tp,0)+i)
        return (NO);
  return (YES);
}

/****************************************************************************/
/*D
   FillRedundantComponentsOfVD - fill the redundant components of a VECDATA_DESC

   SYNOPSIS:
   INT FillRedundantComponentsOfVD (VECDATA_DESC *vd)

   PARAMETERS:
   .  vd - VECDATA_DESC

   DESCRIPTION:
   This function fills the redundant components of a VECDATA_DESC.
   The functions 'ConstructVecOffsets' and 'SetScalVecSettings' are called.

   RETURN VALUE:
   INT
   .n    NUM_OK
   D*/
/****************************************************************************/

INT FillRedundantComponentsOfVD (VECDATA_DESC *vd)
{
  ConstructVecOffsets(VD_NCMPPTR(vd),VD_OFFSETPTR(vd));
  SetCompactTypesOfVec(vd);
  SetScalVecSettings(vd);
  VD_SUCC_COMP(vd) = VDCompsSubsequent(vd);

  return (NUM_OK);
}

/****************************************************************************/
/*D
   GetFirstVector - return first vector for a given multigrid

   SYNOPSIS:
   VECDATA_DESC *GetFirstVector (MULTIGRID *theMG)

   PARAMETERS:
   .  theMG - multigrid

   DESCRIPTION:
   This function returns the first vector for a given multigrid

   RETURN VALUE:
   VECDATA_DESC *
   .n    NULL if there is none
   D*/
/****************************************************************************/

VECDATA_DESC *GetFirstVector (MULTIGRID *theMG)
{
  ENVITEM *item;

  if (ChangeEnvDir("/Multigrids") == NULL) return (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (NULL);
  item = (ENVITEM *)ChangeEnvDir("Vectors");
  if (item == NULL) return (NULL);
  for (item=ENVITEM_DOWN(item); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == VectorVarID)
      return ((VECDATA_DESC *)item);

  return (NULL);
}

/****************************************************************************/
/*D
   GetNextVector - return next vector for a given vector

   SYNOPSIS:
   VECDATA_DESC *GetNextVector (VECDATA_DESC *vd)

   PARAMETERS:
   .  theMG - multigrid

   DESCRIPTION:
   This function returns the next vector for a given vector

   RETURN VALUE:
   VECDATA_DESC *
   .n    NULL if there is no successor
   D*/
/****************************************************************************/

VECDATA_DESC *GetNextVector (VECDATA_DESC *vd)
{
  ENVITEM *item;

  for (item=NEXT_ENVITEM((ENVITEM *)vd); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == VectorVarID)
      return ((VECDATA_DESC *)item);

  return (NULL);
}

static INT GetNewVectorName (MULTIGRID *theMG, char *name)
{
  VECDATA_DESC *vd;
  char buffer[NAMESIZE];
  INT i;

  for (i=0; i<MAX_NAMES; i++) {
    sprintf(buffer,"vec%02d",i);
    for (vd = GetFirstVector(theMG); vd != NULL; vd = GetNextVector(vd))
      if (strcmp(ENVITEM_NAME(vd),buffer) == 0) break;
    if (vd == NULL) break;
  }
  if (i == MAX_NAMES) REP_ERR_RETURN (1);
  strcpy(name,buffer);

  return (NUM_OK);
}

/****************************************************************************/
/*D
   CreateVecDesc - create a vector and fill extra data

   SYNOPSIS:
   VECDATA_DESC *CreateVecDesc (MULTIGRID *theMG, const char *name, const char *compNames,
                                                         const SHORT *NCmpInType, SHORT nId, SHORT *Ident);

   PARAMETERS:
   .  theMG - create vector for this multigrid
   .  name - create vector with this name
   .  NCmpInType - 'VECDATA_DESC' specification
   .  compNames - (optional) vector of component names (in the canonic type order)
               one char each (NULL pointer for no names)
   .  nId - number of comps after identification (maybe

   DESCRIPTION:
   This function creates a 'VECDATA_DESC' and fills its components.

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to 'VECDATA_DESC' environment item
   .n     NULL if failed
   D*/
/****************************************************************************/

VECDATA_DESC *CreateVecDesc (MULTIGRID *theMG, const char *name, const char *compNames,
                             const SHORT *NCmpInType, SHORT nId, SHORT *Ident)
{
  VECDATA_DESC *vd;
  SHORT offset[NVECOFFSETS],*Comp;
  char buffer[NAMESIZE];
  INT i,j,k,tp,ncmp,size;

  if (theMG == NULL)
    REP_ERR_RETURN (NULL);

  if (ChangeEnvDir("/Multigrids") == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir("Vectors") == NULL) {
    MakeEnvItem("Vectors",VectorDirID,sizeof(ENVDIR));
    if (ChangeEnvDir("Vectors") == NULL) REP_ERR_RETURN (NULL);
  }
  if (name != NULL)
    strcpy(buffer,name);
  else if (GetNewVectorName(theMG,buffer)) REP_ERR_RETURN (NULL);
  ConstructVecOffsets(NCmpInType,offset);
  ncmp = offset[NVECTYPES];
  if (ncmp <= 0) REP_ERR_RETURN (NULL);
  size = sizeof(VECDATA_DESC)+(ncmp-1)*sizeof(SHORT);
  vd = (VECDATA_DESC *) MakeEnvItem (buffer,VectorVarID,size);
  if (vd == NULL) REP_ERR_RETURN (NULL);
  if (compNames==NULL)
    memcpy(VM_COMP_NAMEPTR(vd),NoVecNames,MIN(ncmp,MAX_VEC_COMP));
  else
    memcpy(VM_COMP_NAMEPTR(vd),compNames,MIN(ncmp,MAX_VEC_COMP));

  /* fill data in vec data desc */
  VD_MG(vd) = theMG;
  i = 0;
  Comp = VM_COMPPTR(vd);
  for (tp=0; tp<NVECTYPES; tp++) {
    VD_NCMPS_IN_TYPE(vd,tp) = NCmpInType[tp];
    VD_CMPPTR_OF_TYPE(vd,tp) = Comp + offset[tp];
    for (j=0; j<=MAX_NDOF_MOD_32*32-NCmpInType[tp]; j++) {
      if (i >= offset[tp+1]) break;
      if (j*sizeof(DOUBLE) >= FMT_S_VEC_TP(MGFORMAT(theMG),tp))
        REP_ERR_RETURN (NULL);
      if (READ_DR_VEC_FLAG(theMG,tp,j)) continue;
      for (k=1; k<offset[tp+1]-i; k++)
        if (READ_DR_VEC_FLAG(theMG,tp,j+k)) break;
      if ( k<offset[tp+1]-i ) continue;
      Comp[i++] = j;
      SET_DR_VEC_FLAG(theMG,tp,j);
    }
  }
  for (tp=0; tp<NVECOFFSETS; tp++)
    VD_OFFSET(vd,tp) = offset[tp];

  for (tp=0; tp<NVECTYPES; tp++) {
    PRINTDEBUG(np,6,("offset %d comp ",offset[tp]));
    for (i=0; i<VD_NCMPS_IN_TYPE(vd,tp); i++)
      PRINTDEBUG(np,6,(" %d",VD_CMP_OF_TYPE(vd,tp,i)));
  }
  PRINTDEBUG(np,6,("\n"));

  VD_NID(vd) = nId;
  VD_IDENT_PTR(vd) = Ident;

  if (FillRedundantComponentsOfVD(vd))
    REP_ERR_RETURN (NULL);
  VM_LOCKED(vd) = 0;

  return (vd);
}

/****************************************************************************/
/*D
   CreateSubVecDesc - create a sub-vecdesc for a given vector and fill extra data

   SYNOPSIS:
   VECDATA_DESC *CreateSubVecDesc (MULTIGRID *theMG, const char *name,
                                                                const SHORT *NCmpInType, const SHORT *Comps, const char *CompNames)

   PARAMETERS:
   .  theMG - create vector for this multigrid
   .  theVD - given vector
   .  name - create vecdesc with this name (maybe NULL for default name)
   .  NCmpInType - 'VECDATA_DESC' specification
   .  Comps - the Comps of theVD which make up the sub-vecdesc
   .  compNames - (optional) vector of component names (in the canonic type order)
               one char each (NULL pointer for no names)

   DESCRIPTION:
   This function creates a 'VECDATA_DESC' and fills its components.

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to 'VECDATA_DESC' environment item
   .n     NULL if failed
   D*/
/****************************************************************************/

VECDATA_DESC *CreateSubVecDesc (MULTIGRID *theMG, const char *name,
                                const SHORT *NCmpInType, const SHORT *Comps, const char *CompNames)
{
  VECDATA_DESC *vd;
  SHORT offset[NVECOFFSETS];
  INT j,k,tp,ncmp,size;
  char buffer[NAMESIZE];

  if (theMG == NULL)
    REP_ERR_RETURN (NULL);

  if (ChangeEnvDir("/Multigrids") == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir("Vectors") == NULL) REP_ERR_RETURN (NULL);
  ConstructVecOffsets(NCmpInType,offset);
  ncmp = offset[NVECTYPES];
  if (ncmp <= 0) REP_ERR_RETURN (NULL);
  size = sizeof(VECDATA_DESC)+(ncmp-1)*sizeof(SHORT);
  if (name != NULL)
    strcpy(buffer,name);
  else if (GetNewVectorName(theMG,buffer)) REP_ERR_RETURN (NULL);
  vd = (VECDATA_DESC *) MakeEnvItem (buffer,VectorVarID,size);
  if (vd == NULL) REP_ERR_RETURN (NULL);

  /* fill data in vec data desc */
  VD_MG(vd) = theMG;
  strncpy(VM_COMP_NAMEPTR(vd),CompNames,ncmp);
  k = 0;
  for (tp=0; tp<NVECTYPES; tp++) {
    VD_NCMPS_IN_TYPE(vd,tp) = NCmpInType[tp];
    VD_CMPPTR_OF_TYPE(vd,tp) = VM_COMPPTR(vd) + offset[tp];
    for (j=0; j<NCmpInType[tp]; j++) {
      VD_CMP_OF_TYPE(vd,tp,j) = Comps[k++];
    }
  }
  ASSERT(k==offset[NVECTYPES]);
  for (tp=0; tp<NVECOFFSETS; tp++)
    VD_OFFSET(vd,tp) = offset[tp];

  VD_NID(vd) = NO_IDENT;

  if (FillRedundantComponentsOfVD(vd))
    return (NULL);
  VM_LOCKED(vd) = 0;

  return (vd);
}

/****************************************************************************/
/*D
   CombineVecDescs - combines a field of vecdescs to a new vecdesc

   SYNOPSIS:
   VECDATA_DESC *CombineVecDesc (MULTIGRID *theMG, const char *name, const VECDATA_DESC **theVDs,
                                                                const INT nrOfVDs)

   PARAMETERS:
   .  theMG - create vector for this multigrid
   .  name - create vecdesc with this name
   .  theVDs - given vecdescs
   .  nrOfVDs - nrOfVDs vecdescs

   DESCRIPTION:
   This function creates a 'VECDATA_DESC' from a field of vecdescs.
   Components may occur several times! It is up to the user to remove
   redundancies. The resulting vector does lose component names.

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to 'VECDATA_DESC' environment item
   .n     NULL if failed
   D*/
/****************************************************************************/

VECDATA_DESC *CombineVecDesc (MULTIGRID *theMG, const char *name, const VECDATA_DESC **theVDs,
                              const INT nrOfVDs)
{
  VECDATA_DESC *vd;
  SHORT offset;
  INT i,j,k,type,ncmp,size;

  if (theMG == NULL)
    REP_ERR_RETURN (NULL);

  if (ChangeEnvDir("/Multigrids") == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir("Vectors") == NULL) REP_ERR_RETURN (NULL);

  /* compute size of resulting VD */
  ncmp = 0;
  for (i=0; i<nrOfVDs; i++)
  {
    for (type=0; type<NVECTYPES; type++)
      ncmp += VD_NCMPPTR(theVDs[i])[type];
  }
  if (ncmp <= 0) REP_ERR_RETURN (NULL);

  size = sizeof(VECDATA_DESC)+(ncmp-1)*sizeof(SHORT);
  vd = (VECDATA_DESC *) MakeEnvItem (name,VectorVarID,size);
  if (vd == NULL)
    REP_ERR_RETURN (NULL);

  /* fill data in vec data desc */
  strcpy(VM_COMP_NAMEPTR(vd),"");       /* no component names */
  offset=0;
  for (type=0; type<NVECTYPES; type++) {
    VD_OFFSET(vd,type) = offset;
    VD_CMPPTR_OF_TYPE(vd,type) = VM_COMPPTR(vd) + offset;
    k = 0;
    for (i=0; i<nrOfVDs; i++)
    {
      for (j=0; j<VD_NCMPS_IN_TYPE(theVDs[i],type); j++)
        VD_CMP_OF_TYPE(vd,type,k) = VD_CMP_OF_TYPE(theVDs[i],type,j);
    }
    VD_NCMPS_IN_TYPE(vd,type) = k;
    offset += k;
  }
  VD_OFFSET(vd,type) = offset;       /* last one points to the end of the array */

  VD_NID(vd) = NO_IDENT;

  if (FillRedundantComponentsOfVD(vd))
    return (NULL);
  VM_LOCKED(vd) = 0;

  return (vd);
}

/****************************************************************************/
/*D
   VDequal - compare vec data descriptors

   SYNOPSIS:
   INT VDequal (const VECDATA_DESC *vd0, const VECDATA_DESC *vd1)

   PARAMETERS:
   .  vd0 - first  vec data descriptor
   .  vd1 - second vec data descriptor

   DESCRIPTION:
   This function compares two vec data descriptors. It returns YES iff
   they coincide in all components of all types.

   RETURN VALUE:
   INT
   .n      YES if vec data descriptors are equal
   .n      NO  else
 */
/****************************************************************************/

INT VDequal (const VECDATA_DESC *vd0, const VECDATA_DESC *vd1)
{
  INT tp,i,n;
  SHORT *c0,*c1;

  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_NCMPS_IN_TYPE(vd0,tp)==VD_NCMPS_IN_TYPE(vd1,tp))
    {
      n  = VD_NCMPS_IN_TYPE(vd0,tp);
      if (n<=0)
        continue;

      c0 = VD_CMPPTR_OF_TYPE(vd0,tp);
      c1 = VD_CMPPTR_OF_TYPE(vd1,tp);
      for (i=0; i<n; i++)
        if (c0[i]!=c1[i])
          return (NO);
    }
    else
      return (NO);

  return (YES);
}

/****************************************************************************/
/*D
   AllocVDfromVD - dynamic vector allocation

   SYNOPSIS:
   INT AllocVDfromVD (MULTIGRID *theMG, INT fl, INT tl,
   VECDATA_DESC *template_desc, VECDATA_DESC **new_desc);

   PARAMETERS:
   .  theMG - create vector for this multigrid
   .  fl - from level
   .  tl - to level
   .  template_desc - template vector
   .  new_desc - new vector

   DESCRIPTION:
   This function allocates a new vector.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

static INT CompVecDesc (const VECDATA_DESC *vd, const SHORT *NCmpsInType)
{
  INT tp;

  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_NCMPS_IN_TYPE(vd,tp) != NCmpsInType[tp])
      return(1);

  return(0);
}

static INT IsVecDescAlloc (MULTIGRID *theMG, INT fl, INT tl, const VECDATA_DESC *vd)
{
  GRID *theGrid;
  INT i,j,tp;

  if (vd == NULL) return(NO);

  /* are the components of vd free */
  for (i=fl; i<=tl; i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (tp=0; tp<NVECTYPES; tp++)
      for (j=0; j<VD_NCMPS_IN_TYPE(vd,tp); j++)
        if (READ_DR_VEC_FLAG(theGrid,tp,VD_CMP_OF_TYPE(vd,tp,j))==NO)
          return(NO);
  }

  return (YES);
}

static INT GetVecDescAllocLevels (MULTIGRID *theMG, const VECDATA_DESC *vd, INT AllocLev[MAXLEVEL])
{
  INT lev;

  for (lev=0; lev<MAXLEVEL; lev++) AllocLev[lev] = NO;

  if (vd == NULL) REP_ERR_RETURN(1);

  for (lev=0; lev<=TOPLEVEL(theMG); lev++)
    AllocLev[lev] = IsVecDescAlloc(theMG,lev,lev,vd);

  return (NUM_OK);
}

static INT AllocVecDesc (MULTIGRID *theMG, INT fl, INT tl, const VECDATA_DESC *vd)
{
  GRID *theGrid;
  INT i,j,tp;

  if (vd == NULL) return(1);            /* cannot allocate */

  /* are the components of vd free */
  for (i=fl; i<=tl; i++) {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (tp=0; tp<NVECTYPES; tp++)
      for (j=0; j<VD_NCMPS_IN_TYPE(vd,tp); j++)
        if (READ_DR_VEC_FLAG(theGrid,tp,VD_CMP_OF_TYPE(vd,tp,j)))
          return(1);                                    /* NO */
  }

  /* allocate vd without changing components */
  for (i=fl; i<=tl; i++) {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (tp=0; tp<NVECTYPES; tp++)
      for (j=0; j<VD_NCMPS_IN_TYPE(vd,tp); j++)
        SET_DR_VEC_FLAG(theGrid,tp,VD_CMP_OF_TYPE(vd,tp,j));
  }
  for (tp=0; tp<NVECTYPES; tp++)
    for (j=0; j<VD_NCMPS_IN_TYPE(vd,tp); j++)
      SET_DR_VEC_FLAG(theMG,tp,VD_CMP_OF_TYPE(vd,tp,j));

  return (NUM_OK);
}

/****************************************************************************/
/*D
   AllocVDfromNCmp - dynamically allocate vector descriptor from given NCmpInType and compNames

   SYNOPSIS:
   INT AllocVDfromNCmp (MULTIGRID *theMG, INT fl, INT tl,
                                   const SHORT *NCmpInType, const char *compNames, VECDATA_DESC **new_desc)

   PARAMETERS:
   .  theMG - allocate descriptor for this multigrid
   .  fl - from level
   .  tl - to level
   .  NCmpInType - rows per type
   .  compNames  - names of comps (may be NULL)
   .  new_desc   - handle for descriptor

   DESCRIPTION:
   This function allocates a vector descriptor from
   given NCmpInType and compNames. If there is no existing
   decriptor a new one is created.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT AllocVDfromNCmp (MULTIGRID *theMG, INT fl, INT tl,
                     const SHORT *NCmpInType, const char *compNames, VECDATA_DESC **new_desc)
{
  VECDATA_DESC *vd;

  if (*new_desc != NULL)
    if (VM_LOCKED(*new_desc))
      return (NUM_OK);
  if (AllocVecDesc(theMG,fl,tl,*new_desc)) {
    for (vd = GetFirstVector(theMG); vd != NULL; vd = GetNextVector(vd)) {
      if (VM_LOCKED(vd)) continue;
      if (CompVecDesc(vd,NCmpInType)) continue;
      if (!AllocVecDesc(theMG,fl,tl,vd)) {
        *new_desc = vd;
        return (NUM_OK);
      }
    }
    *new_desc = CreateVecDesc(theMG,NULL,compNames,
                              NCmpInType,NO_IDENT,NULL);
    if (*new_desc == NULL) {
      PrintErrorMessage('E',"AllocVDfromNCmp","cannot create VecDesc\n");
      REP_ERR_RETURN (1);
    }
    if (AllocVecDesc(theMG,fl,tl,*new_desc)) {
      PrintErrorMessage('E',"AllocVDfromNCmp","cannot allocate VecDesc\n");
      REP_ERR_RETURN(1);
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   AllocVDFromVD - dynamically allocate vector descriptor from given vector descriptor

   SYNOPSIS:
   INT AllocVDFromVD (MULTIGRID *theMG, INT fl, INT tl,
                                   const VECDATA_DESC *vd, VECDATA_DESC **new_desc)

   PARAMETERS:
   .  theMG - allocate descriptor for this multigrid
   .  fl - from level
   .  tl - to level
   .  vd - given vector descriptor
   .  new_desc   - handle for descriptor

   DESCRIPTION:
   This function allocates a vector descriptor from a
   given vector descriptors.
   If there is no existing decriptor a new one is created.

   SEE ALSO:
   AllocVDfromNCmp

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT AllocVDFromVD (MULTIGRID *theMG, INT fl, INT tl,
                   const VECDATA_DESC *vd, VECDATA_DESC **new_desc)
{
  if (AllocVDfromNCmp(theMG,fl,tl,vd->NCmpInType,vd->compNames,new_desc))
    REP_ERR_RETURN(1);

  VD_NID(*new_desc) = VD_NID(vd);
  VD_IDENT_PTR(*new_desc) = VD_IDENT_PTR(vd);

  return (0);
}

/****************************************************************************/
/*D
   LockVD - protect vector against removal or deallocation

   SYNOPSIS:
   INT LockVD (VECDATA_DESC *vd)

   PARAMETERS:
   .  vd - vector descriptor

   DESCRIPTION:
   This function locks a vector against removal or deallocation.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT LockVD (VECDATA_DESC *vd)
{
  VM_LOCKED(vd) = VM_IS_LOCKED;
  return (0);
}

/****************************************************************************/
/*D
   TransmitLockStatusVD - ...

   SYNOPSIS:
   INT TransmitLockStatusVD (const VECDATA_DESC *vd, VECDATA_DESC *svd)

   PARAMETERS:
   .  vd  - vector descriptor
   .  svd - sub vector descriptor

   DESCRIPTION:
   This function ...

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT TransmitLockStatusVD (const VECDATA_DESC *vd, VECDATA_DESC *svd)
{
  if (!VM_LOCKED(vd) && VM_LOCKED(svd))
    REP_ERR_RETURN(1);
  VM_LOCKED(svd) = VM_LOCKED(vd);

  return (0);
}

/****************************************************************************/
/*D
   FreeVD - dynamic vector deallocation

   SYNOPSIS:
   INT FreeVD (MULTIGRID *theMG, INT fl, INT tl, VECDATA_DESC *vd);

   PARAMETERS:
   .  theMG - create vector for this multigrid
   .  fl - from level
   .  tl - to level
   .  vd - vector descriptor

   DESCRIPTION:
   This function deallocates a vector.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT FreeVD (MULTIGRID *theMG, INT fl, INT tl, VECDATA_DESC *vd)
{
  GRID *theGrid;
  INT i,j,tp;

  if (vd==NULL) return (NUM_OK);
  dset(theMG,fl,tl,ALL_VECTORS,vd,0.0);
  if (VM_LOCKED(vd)) return (NUM_OK);
  PRINTDEBUG(np,2,(" FreeVD %s from %d to %d\n",
                   ENVITEM_NAME(vd),fl,tl));
  for (i=fl; i<=tl; i++) {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (tp=0; tp<NVECTYPES; tp++)
      for (j=0; j<VD_NCMPS_IN_TYPE(vd,tp); j++)
        CLEAR_DR_VEC_FLAG(theGrid,tp,VD_CMP_OF_TYPE(vd,tp,j));
  }

  for (i=BOTTOMLEVEL(theMG); i<=TOPLEVEL(theMG); i++) {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (tp=0; tp<NVECTYPES; tp++)
      for (j=0; j<VD_NCMPS_IN_TYPE(vd,tp); j++)
        if (READ_DR_VEC_FLAG(theGrid,tp,VD_CMP_OF_TYPE(vd,tp,j)))
          return(NUM_OK);
  }
  for (tp=0; tp<NVECTYPES; tp++)
    for (j=0; j<VD_NCMPS_IN_TYPE(vd,tp); j++)
      CLEAR_DR_VEC_FLAG(theMG,tp,VD_CMP_OF_TYPE(vd,tp,j));


  return (NUM_OK);
}

/****************************************************************************/
/*D
   InterpolateVDAllocation - dynamic vector allocation on new level

   SYNOPSIS:
   INT InterpolateVDAllocation (MULTIGRID *theMG, VECDATA_DESC *vd);

   PARAMETERS:
   .  theMG -  multigrid
   .  vd - vector descriptor

   DESCRIPTION:
   This function allocates a vector on a new level.

   RETURN VALUE:
   INT
   .n      NUM_OK if ok
   .n      NUM_ERROR if error occurred
 */
/****************************************************************************/

INT InterpolateVDAllocation (MULTIGRID *theMG, VECDATA_DESC *vd)
{
  GRID *theGrid;
  INT j,tp,tl;

  if (vd==NULL) return (NUM_OK);
  if (VM_LOCKED(vd)) return (NUM_OK);
  tl = TOPLEVEL(theMG);
  if (tl < 1) return (NUM_OK);

  PRINTDEBUG(np,2,(" InterpolateVDAllocation %s\n",ENVITEM_NAME(vd)));

  theGrid = GRID_ON_LEVEL(theMG,tl);
  for (tp=0; tp<NVECTYPES; tp++)
    for (j=0; j<VD_NCMPS_IN_TYPE(vd,tp); j++) {
      if (READ_DR_VEC_FLAG(theGrid,tp,VD_CMP_OF_TYPE(vd,tp,j)))
        return(NUM_ERROR);
      SET_DR_VEC_FLAG(theGrid,tp,VD_CMP_OF_TYPE(vd,tp,j));
    }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   DisposeVD - remove VECDATA_DESC from objects of multigrid

   SYNOPSIS:
   INT DisposeVD (VECDATA_DESC *vd)

   PARAMETERS:
   .  vd - vector descriptor

   DESCRIPTION:
   This function removes a vector descriptor from the objects of the multigrid. The part
   of environment memory it occupies is freed. The functions calls
   'RemoveEnvItem'.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT DisposeVD (VECDATA_DESC *vd)
{
  if (vd==NULL) REP_ERR_RETURN (NUM_ERROR);
  if (VM_LOCKED(vd)) REP_ERR_RETURN (NUM_ERROR);

  ENVITEM_LOCKED(vd) = 0;

  if (ChangeEnvDir("/Multigrids") == NULL) REP_ERR_RETURN (0);
  if (ChangeEnvDir(ENVITEM_NAME(VD_MG(vd))) == NULL) REP_ERR_RETURN (0);
  if (ChangeEnvDir("Vectors") == NULL) REP_ERR_RETURN (0);

  if (RemoveEnvItem((ENVITEM*)vd))
    REP_ERR_RETURN (0);

  return (0);
}

/****************************************************************************/
/*D
   DisplayVecDataDesc - Display VECDATA_DESC entries

   SYNOPSIS:
   INT DisplayVecDataDesc (const VECDATA_DESC *vd, char *buffer, INT modifiers)

   PARAMETERS:
   .  fmt - associated format for names of abstract types
   .  vd - VECDATA_DESC to display
   .  buffer - print here
   .  modifiers - modifier flags

   DESCRIPTION:
   This function displays the entries of a VECDATA_DESC: comp-names, comp-positions etc.

   RETURN VALUE:
   INT
   .n      0: ok
   .n      else: error
   D*/
/****************************************************************************/

static INT lev2str (const levels[MAXLEVEL], char *list)
{
  INT i,f,t,p;

  p = 0;
  for (i=0; i<MAXLEVEL; i++)
  {
    /* skip NOs */
    while (!levels[i] && i<MAXLEVEL) i++;

    if (i>=MAXLEVEL)
      if (p==0)
        return (1);
      else
        break;

    f = i;

    /* skip NOs */
    while (levels[i] && i<MAXLEVEL) i++;

    t = i-1;

    switch (t-f)
    {
    case 0 :
      p += sprintf(list+p,"%d,",f);
      break;
    case 1 :
      p += sprintf(list+p,"%d,%d,",f,t);
      break;
    default :
      p += sprintf(list+p,"%d-%d,",f,t);
    }
  }
  list[p-1] = '\0';

  return (0);
}

INT DisplayVecDataDesc (const VECDATA_DESC *vd, INT modifiers, char *buffer)
{
  const FORMAT *fmt;
  const SHORT *offset;
  const char *cn;
  INT rt,i;

  if (vd==NULL) REP_ERR_RETURN (1);

  buffer += sprintf(buffer,"vector data descriptor '%s'\n",ENVITEM_NAME(vd));

  fmt = MGFORMAT(VD_MG(vd));
  cn = VM_COMP_NAMEPTR(vd);
  offset = VD_OFFSETPTR(vd);
  for (rt=0; rt<NVECTYPES; rt++)
    if (VD_ISDEF_IN_TYPE(vd,rt))
    {
      buffer += sprintf(buffer,"-------\n");
      for (i=0; i<VD_NCMPS_IN_TYPE(vd,rt); i++)
        buffer += sprintf(buffer,"%c  %c %2d\n",(i) ? ' ' : FMT_VTYPE_NAME(fmt,rt),cn[offset[rt]+i],VD_CMP_OF_TYPE(vd,rt,i));
    }
  buffer += sprintf(buffer,"-------\n");

  if (READ_FLAG(modifiers,SCAL_PROP))
    if (VD_IS_SCALAR(vd))
    {
      buffer += sprintf(buffer,"\ndescriptor is scalar:\n");
      buffer += sprintf(buffer,"  comp %2d\n",VD_SCALCMP(vd));
      buffer += sprintf(buffer,"  mask %2d\n",VD_SCALTYPEMASK(vd));
    }

  if (READ_FLAG(modifiers,ALLOC_STAT))
  {
    if (VM_LOCKED(vd))
      buffer += sprintf(buffer,"descriptor is locked\n");
    else
    {
      INT levels[MAXLEVEL];
      char LevelList[MAXLEVEL];

      if (GetVecDescAllocLevels(VD_MG(vd),vd,levels))
        REP_ERR_RETURN (1);
      if (lev2str(levels,LevelList))
        buffer += sprintf(buffer,"descriptor is not allocated\n");
      else
        buffer += sprintf(buffer,"descriptor is allocated on levels [%s]\n",LevelList);
    }
  }
  buffer += sprintf(buffer,"\n");

  return (NUM_OK);
}

/****************************************************************************/
/*D
   GetVecDataDescByName - find vector data desciptor

   SYNOPSIS:
   VECDATA_DESC *GetVecDataDescByName (MULTIGRID *theMG, char *name);

   PARAMETERS:
   .  theMG - create vector for this multigrid
   .  name - name of a vector

   DESCRIPTION:
   This function finds a vector by name.

   RETURN VALUE:
   VECDATA_DESC *
   .n      pointer to the vector
   .n      NULL if there is no vector of this name in the multigrid
   D*/
/****************************************************************************/

VECDATA_DESC *GetVecDataDescByName (const MULTIGRID *theMG, char *name)
{
  if (ChangeEnvDir("/Multigrids") == NULL) return (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (NULL);
  return((VECDATA_DESC *) SearchEnv(name,"Vectors",
                                    VectorVarID,VectorDirID));
}

/****************************************************************************/
/*D
        VDinterfaceDesc - an interface VECDATA_DESC is created

        SYNOPSIS:
        INT VDinterfaceDesc (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi)

    PARAMETERS:
   .   vd			- make a sub desc of this VECDATA_DESC
   .   vds			- an existing sub desc of vd
   .   vdi			- handle to new interface desc

        DESCRIPTION:
        This function creates a sub descriptor to a given VECDATA_DESC such that all components
        from vds are taken of types in which vds is defined but where vds is a true subset of vd.

        RETURN VALUE:
        INT
   .n   0: ok
   .n      n: if an error occured
   D*/
/****************************************************************************/

INT VDinterfaceDesc (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi)
{
  SHORT SubComp[MAX_VEC_COMP],SubNCmp[NVECTYPES];
  INT i,k,n,ns,tp;
  char SubName[MAX_VEC_COMP],buffer[NAMESIZE];

  /* generate name and see if desc already exists */
  strcpy(buffer,ENVITEM_NAME(vds));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,"i");
  *vdi = GetVecDataDescByName(MD_MG(vd),buffer);
  if (*vdi != NULL) {
    if (TransmitLockStatusVD(vds,*vdi))
      REP_ERR_RETURN(1);
    return(0);
  }

  k = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vds,tp))
    {
      if (!VD_ISDEF_IN_TYPE(vd,tp))
        REP_ERR_RETURN (1);

      n  = VD_NCMPS_IN_TYPE(vd,tp);
      ns = VD_NCMPS_IN_TYPE(vds,tp);
      if (ns<n)
      {
        /* copy all components from vds */
        for (i=0; i<ns; i++)
        {
          ASSERT(k<MAX_MAT_COMP);

          SubComp[k] = VD_CMP_OF_TYPE(vds,tp,i);
          SubName[k] = VM_COMP_NAME(vds,VD_OFFSET(vds,tp)+i);
          k++;
        }
        SubNCmp[tp] = ns;
      }
      else if (ns==n)
      {
        /* no components here */
        SubNCmp[tp] = 0;
      }
      else
        /* vd does not contain vds */
        REP_ERR_RETURN (1);
    }
    else
      /* no components here */
      SubNCmp[tp] = 0;

  *vdi = CreateSubVecDesc(VD_MG(vd),buffer,SubNCmp,SubComp,SubName);
  if (*vdi == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusVD(vd,*vdi))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/*D
        VDinterfaceCoDesc - an interface co-VECDATA_DESC is created

        SYNOPSIS:
        INT VDinterfaceCoDesc (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi)

    PARAMETERS:
   .   vd			- make a sub desc of this VECDATA_DESC
   .   vds			- an existing sub desc of vd
   .   vdi			- handle to new interface desc

        DESCRIPTION:
        This function creates a sub descriptor to a given VECDATA_DESC such that all components
        from vd are taken of types in which vds is defined but where vds is a true subset of vd.
        Only components of vd are taken that are not in vds.

        RETURN VALUE:
        INT
   .n   0: ok
   .n      n: if an error occured
   D*/
/****************************************************************************/

INT VDinterfaceCoDesc (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi)
{
  SHORT SubComp[MAX_VEC_COMP],SubNCmp[NVECTYPES];
  INT i,j,k,n,ns,tp,cmp,ncmp;
  char SubName[MAX_VEC_COMP],buffer[NAMESIZE];

  /* generate name and see if desc already exists */
  strcpy(buffer,ENVITEM_NAME(vds));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,"ico");
  *vdi = GetVecDataDescByName(MD_MG(vd),buffer);
  if (*vdi != NULL) {
    if (TransmitLockStatusVD(vds,*vdi))
      REP_ERR_RETURN(1);
    return(0);
  }

  k = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vds,tp))
    {
      if (!VD_ISDEF_IN_TYPE(vd,tp))
        REP_ERR_RETURN (1);

      n  = VD_NCMPS_IN_TYPE(vd,tp);
      ns = VD_NCMPS_IN_TYPE(vds,tp);
      if (ns<n)
      {
        /* copy all components from vd not in vds */
        ncmp = 0;
        for (i=0; i<n; i++)
        {
          cmp = VD_CMP_OF_TYPE(vd,tp,i);
          for (j=0; j<ns; j++)
            if (VD_CMP_OF_TYPE(vds,tp,j)==cmp)
              break;
          if (j<ns)
            /* cmp contained in vds */
            continue;

          ASSERT(k<MAX_MAT_COMP);

          SubComp[k] = cmp;
          SubName[k] = VM_COMP_NAME(vd,VD_OFFSET(vd,tp)+i);
          k++;
          ncmp++;
        }
        SubNCmp[tp] = ncmp;
      }
      else if (ns==n)
      {
        /* no components here */
        SubNCmp[tp] = 0;
      }
      else
        /* vd does not contain vds */
        REP_ERR_RETURN (1);
    }
    else
      /* no components here */
      SubNCmp[tp] = 0;

  *vdi = CreateSubVecDesc(VD_MG(vd),buffer,SubNCmp,SubComp,SubName);
  if (*vdi == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusVD(vd,*vdi))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/*D
        VDCoDesc - a co-VECDATA_DESC is created

        SYNOPSIS:
        INT VDCoDesc (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi)

    PARAMETERS:
   .   vd			- make a sub desc of this VECDATA_DESC
   .   vds			- an existing sub desc of vd
   .   vdi			- handle to new interface desc

        DESCRIPTION:
        This function creates a sub descriptor to a given VECDATA_DESC vd such that it is
        the complement of another given descriptor vds (which is also sub to vd).

        RETURN VALUE:
        INT
   .n   0: ok
   .n      n: if an error occured
   D*/
/****************************************************************************/

INT VDCoDesc (const VECDATA_DESC *vd, const VECDATA_DESC *vds, VECDATA_DESC **vdi)
{
  SHORT SubComp[MAX_VEC_COMP],SubNCmp[NVECTYPES];
  INT i,j,k,n,ns,tp,cmp,ncmp;
  char SubName[MAX_VEC_COMP],buffer[NAMESIZE];

  /* generate name and see if desc already exists */
  strcpy(buffer,ENVITEM_NAME(vds));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,"co");
  *vdi = GetVecDataDescByName(MD_MG(vd),buffer);
  if (*vdi != NULL) {
    if (TransmitLockStatusVD(vds,*vdi))
      REP_ERR_RETURN(1);
    return(0);
  }

  k = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
    {
      if (VD_ISDEF_IN_TYPE(vds,tp))
      {

        n  = VD_NCMPS_IN_TYPE(vd,tp);
        ns = VD_NCMPS_IN_TYPE(vds,tp);
        if (ns<n)
        {
          /* copy all components from vd not in vds */
          ncmp = 0;
          for (i=0; i<n; i++)
          {
            cmp = VD_CMP_OF_TYPE(vd,tp,i);
            for (j=0; j<ns; j++)
              if (VD_CMP_OF_TYPE(vds,tp,j)==cmp)
                break;
            if (j<ns)
              /* cmp contained in vds */
              continue;

            ASSERT(k<MAX_MAT_COMP);

            SubComp[k] = cmp;
            SubName[k] = VM_COMP_NAME(vd,VD_OFFSET(vd,tp)+i);
            k++;
            ncmp++;
          }
          SubNCmp[tp] = ncmp;
        }
        else if (ns==n)
        {
          /* no components here */
          SubNCmp[tp] = 0;
        }
        else
          /* vd does not contain vds */
          REP_ERR_RETURN (1);
      }
      else
      {
        /* copy all components from vd not in vds */
        n  = VD_NCMPS_IN_TYPE(vd,tp);
        for (i=0; i<n; i++)
        {
          ASSERT(k<MAX_MAT_COMP);

          SubComp[k] = VD_CMP_OF_TYPE(vd,tp,i);
          SubName[k] = VM_COMP_NAME(vd,VD_OFFSET(vd,tp)+i);
          k++;
        }
        SubNCmp[tp] = n;
      }
    }
    else
      /* no components here */
      SubNCmp[tp] = 0;

  *vdi = CreateSubVecDesc(VD_MG(vd),buffer,SubNCmp,SubComp,SubName);
  if (*vdi == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusVD(vd,*vdi))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/*D
   VD_ncmps_in_otype_mod - return number of comps in object if unique, -1 else

   SYNOPSIS:
   INT VD_ncmps_in_otype_mod (const VECDATA_DESC *vd, INT otype, INT mode)

   PARAMETERS:
   .  vd - data decsriptor
   .  otype - object type
   .  mode - STRICT or NON_STRICT

   DESCRIPTION:
   This function checks whether the number of components described in 'otype'
   is the same for all vtypes using objects of 'otype' and returns it.
   If the number is not unique a -1 is returned.
   If not the whole domain is covered, a -2 is returned.
   The uniqueness of comps is not checked here. If mode is STRICT vectors are
   required in ALL vobjects of otype.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'otype'
   are not defined at all!

   RETURN VALUE:
   INT
   .n      number of components in objects of 'otype'
   .n      -1 if not unique
   .n      -2 if not the whole domain is covered
   D*/
/****************************************************************************/

INT VD_ncmps_in_otype_mod (const VECDATA_DESC *vd, INT otype, INT mode)
{
  FORMAT *fmt;
  INT tp,otp,ncmp,parts,i,n;

  fmt = MGFORMAT(VD_MG(vd));
  otp = 1<<otype;
  ncmp = parts = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
      if (otp & FMT_T2O(fmt,tp))
      {
        if (ncmp==0)
          ncmp = VD_NCMPS_IN_TYPE(vd,tp);
        else
        if (VD_NCMPS_IN_TYPE(vd,tp)!=ncmp)
          REP_ERR_RETURN (-1);
        parts |= FMT_T2P(fmt,tp);
      }

  if (mode==STRICT)
  {
    /* check whether all parts are covered */
    n = BVPD_NPARTS(MG_BVPD(VD_MG(vd)));
    for (i=0; i<n; i++)
      if (!(parts & (1<<i)))
        REP_ERR_RETURN (-2);
  }
  else if (mode!=NON_STRICT)
    REP_ERR_RETURN (-3);

  return (ncmp);
}

/****************************************************************************/
/*D
   VD_cmp_of_otype_mod - return comp in object if unique, -1 else

   SYNOPSIS:
   INT VD_cmp_of_otype_mod (const VECDATA_DESC *vd, INT otype, INT i, INT mode)

   PARAMETERS:
   .  vd - data decsriptor
   .  otype - object type
   .  mode - STRICT or NON_STRICT

   DESCRIPTION:
   This function checks whether the offset of component i described in 'otype'
   is the same for all vtypes using objects of 'otype' and returns it.
   If the offset is not unique a -1 is returned. If mode is STRICT vectors are
   required in ALL vobjects of otype.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'otype'
   are not defined at all!

   RETURN VALUE:
   INT
   .n      number of components in objects of 'otype'
   .n      -1 if not unique
   D*/
/****************************************************************************/

INT VD_cmp_of_otype_mod (const VECDATA_DESC *vd, INT otype, INT i, INT mode)
{
  FORMAT *fmt;
  INT tp,otp,ncmp,off,parts,j,n;

  fmt = MGFORMAT(VD_MG(vd));
  otp = 1<<otype;
  ncmp = off = parts = 0;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
      if (otp & FMT_T2O(fmt,tp))
      {
        if (ncmp==0)
        {
          ncmp = VD_NCMPS_IN_TYPE(vd,tp);
          off  = VD_CMP_OF_TYPE(vd,tp,i);
          if (i>=ncmp)
            REP_ERR_RETURN (-1);
        }
        else
        {
          if (VD_NCMPS_IN_TYPE(vd,tp)!=ncmp)
            REP_ERR_RETURN (-1);
          if (VD_CMP_OF_TYPE(vd,tp,i)!=off)
            REP_ERR_RETURN (-1);
        }
        parts |= FMT_T2P(fmt,tp);
      }

  if (mode==STRICT)
  {
    /* check whether all parts are covered */
    n = BVPD_NPARTS(MG_BVPD(VD_MG(vd)));
    for (j=0; j<n; j++)
      if (!(parts & (1<<j)))
        REP_ERR_RETURN (-2);
  }
  else if (mode!=NON_STRICT)
    REP_ERR_RETURN (-3);

  return (off);
}

/****************************************************************************/
/*D
   VD_ncmp_cmpptr_of_otype_mod - return comp in object if unique, -1 else

   SYNOPSIS:
   SHORT *VD_ncmp_cmpptr_of_otype_mod (const VECDATA_DESC *vd, INT otype, INT *ncomp, INT mode)

   PARAMETERS:
   .  vd - data decsriptor
   .  otype - object type
   .  ncomp - number of components (may be NULL)
   .  mode - STRICT or NON_STRICT

   DESCRIPTION:
   This function checks whether all components described in 'otype'
   are the same for all vtypes using objects of 'otype' and returns a component pointer.
   If the components are not unique a NULL is returned. If mode is STRICT vectors are
   required in ALL vobjects of otype.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'otype'
   are not defined at all!

   RETURN VALUE:
   SHORT *
   .n      number of components in objects of 'otype'
   .n      NULL if not unique
   D*/
/****************************************************************************/

SHORT *VD_ncmp_cmpptr_of_otype_mod (const VECDATA_DESC *vd, INT otype, INT *ncomp, INT mode)
{
  FORMAT *fmt;
  SHORT *cptr;
  INT tp,otp,ncmp,i,parts,n;

  if (ncomp!=NULL) *ncomp = -1;

  fmt = MGFORMAT(VD_MG(vd));
  otp = 1<<otype;
  ncmp = parts = 0;
  cptr = NULL;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
      if (otp & FMT_T2O(fmt,tp))
      {
        if (ncmp==0)
        {
          ncmp = VD_NCMPS_IN_TYPE(vd,tp);
          cptr = VD_CMPPTR_OF_TYPE(vd,tp);
        }
        else
        {
          if (VD_NCMPS_IN_TYPE(vd,tp)!=ncmp)
            REP_ERR_RETURN (NULL);
          for (i=0; i<ncmp; i++)
            if (VD_CMP_OF_TYPE(vd,tp,i)!=cptr[i])
              REP_ERR_RETURN (NULL);
        }
        parts |= FMT_T2P(fmt,tp);
      }

  if (mode==STRICT)
  {
    /* check whether all parts are covered */
    n = BVPD_NPARTS(MG_BVPD(VD_MG(vd)));
    for (i=0; i<n; i++)
      if (!(parts & (1<<i)))
        REP_ERR_RETURN (NULL);
  }
  else if (mode!=NON_STRICT)
    REP_ERR_RETURN (NULL);

  if (ncomp!=NULL) *ncomp = ncmp;

  return (cptr);
}

/****************************************************************************/
/*D
   VDusesVOTypeOnly - check whether only one vector object type is used

   SYNOPSIS:
   INT VDusesVOTypeOnly (const VECDATA_DESC *vd, INT votype)

   PARAMETERS:
   .  vd - data decsriptor
   .  votype - vector object type

   DESCRIPTION:
   This function checks whether only one vector object type is used by the
   VECDATA_DESC.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'votype'
   are not defined at all!

   RETURN VALUE:
   SHORT *
   .n      number of components in objects of 'otype'
   .n      NULL if not unique
   D*/
/****************************************************************************/

INT VDusesVOTypeOnly (const VECDATA_DESC *vd, INT votype)
{
  FORMAT *fmt;
  INT tp,otp;

  fmt = MGFORMAT(VD_MG(vd));
  otp = 1<<votype;
  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_ISDEF_IN_TYPE(vd,tp))
      if (otp!=FMT_T2O(fmt,tp))
        return(NO);
  return (YES);
}

/****************************************************************************/
/****************************************************************************/
/*			here follows matrix stuff										*/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*D
   ConstructMatOffsets - Calculate offsets for MAT_SCALARs

   SYNOPSIS:
   INT ConstructMatOffsets (const SHORT *RowsInType, const SHORT *ColsInType, SHORT *offset)

   PARAMETERS:
   .  RowsInType - the numbers of row components is used for the calculation of offsets
   .  ColsInType - the numbers of col components is used for the calculation of offsets
   .  offset - array of length NMATOFFSETS (!)

   DESCRIPTION:
   This function calculates offsets for MAT_SCALARs.

   RETURN VALUE:
   INT
   .n    NUM_OK
   D*/
/****************************************************************************/

INT ConstructMatOffsets (const SHORT *RowsInType, const SHORT *ColsInType, SHORT *offset)
{
  INT type;

  offset[0] = 0;
  for (type=0; type<NMATTYPES; type++)
    offset[type+1] = offset[type] + RowsInType[type]*ColsInType[type];

  return (NUM_OK);
}

/****************************************************************************/
/*D
   ConstructMatOffsetsAlt - Calculate offsets for MAT_SCALARs

   SYNOPSIS:
   INT ConstructMatOffsetsAlt (const SHORT *CmpsInType, SHORT *offset)


   PARAMETERS:
   .  CmpsInType - the numbers of components is used for the calculation of offsets
   .  offset - array of length NMATOFFSETS (!)

   DESCRIPTION:
   This function calculates offsets for MAT_SCALARs.

   RETURN VALUE:
   INT
   .n    NUM_OK
   D*/
/****************************************************************************/

INT ConstructMatOffsetsAlt (const SHORT *CmpsInType, SHORT *offset)
{
  INT type;

  offset[0] = 0;
  for (type=1; type<NMATOFFSETS; type++)
    offset[type] = offset[type-1] + CmpsInType[type-1];

  return (NUM_OK);
}

/****************************************************************************/
/*D
   SetScalMatSettings - fill the scalar settings components of a MATDATA_DESC

   SYNOPSIS:
   INT SetScalMatSettings (MATDATA_DESC *md)

   PARAMETERS:
   .  md - fill this scalar settings components

   DESCRIPTION:
   This function fills the scalar settings components of a MATDATA_DESC.

   RETURN VALUE:
   INT
   .n    NUM_OK
   D*/
/****************************************************************************/

static INT SetScalMatSettings (MATDATA_DESC *md)
{
  INT mtp;

  MD_IS_SCALAR(md) = FALSE;

  /* check number of components per type */
  for (mtp=0; mtp<NMATTYPES; mtp++)
    if (MD_ISDEF_IN_MTYPE(md,mtp))
    {
      if ((MD_ROWS_IN_MTYPE(md,mtp)!=1) || (MD_COLS_IN_MTYPE(md,mtp)!=1))
        return (NUM_OK);                                                        /* no scalar */
      else
        MD_SCALCMP(md) = MD_MCMP_OF_MTYPE(md,mtp,0);
    }

  /* check location of components per type */
  MD_SCAL_RTYPEMASK(md) = MD_SCAL_CTYPEMASK(md) = 0;
  for (mtp=0; mtp<NMATTYPES; mtp++)
    if (MD_ISDEF_IN_MTYPE(md,mtp))
    {
      MD_SCAL_RTYPEMASK(md) |= 1<<MTYPE_RT(mtp);
      MD_SCAL_CTYPEMASK(md) |= 1<<MTYPE_CT(mtp);
      if (MD_SCALCMP(md)!=MD_MCMP_OF_MTYPE(md,mtp,0))
        return (NUM_OK);                                                        /* no scalar */
    }

  MD_IS_SCALAR(md) = TRUE;

  return (NUM_OK);
}

static INT SetCompactTypesOfMat (MATDATA_DESC *md)
{
  FORMAT *fmt;
  INT rt,ct;

  /* fill bitwise fields */
  fmt = MGFORMAT(MD_MG(md));
  MD_ROW_DATA_TYPES(md) = MD_COL_DATA_TYPES(md) =
                            MD_ROW_OBJ_USED(md) = MD_COL_OBJ_USED(md) = 0;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
      {
        MD_ROW_DATA_TYPES(md) |= BITWISE_TYPE(rt);
        MD_COL_DATA_TYPES(md) |= BITWISE_TYPE(ct);
        MD_ROW_OBJ_USED(md)   |= FMT_T2O(fmt,rt);
        MD_COL_OBJ_USED(md)   |= FMT_T2O(fmt,ct);
      }
  return (0);
}

static INT MDCompsSubsequent (const MATDATA_DESC *md)
{
  INT tp,i;

  for (tp=0; tp<NMATTYPES; tp++)
    for (i=0; i<MD_NCMPS_IN_MTYPE(md,tp); i++)
      if (MD_MCMP_OF_MTYPE(md,tp,i)!=MD_MCMP_OF_MTYPE(md,tp,0)+i)
        return (NO);
  return (YES);
}

/****************************************************************************/
/*D
   FillRedundantComponentsOfMD - fill the redundant components of a MATDATA_DESC

   SYNOPSIS:
   INT FillRedundantComponentsOfMD (MATDATA_DESC *md)

   PARAMETERS:
   .  md - MATDATA_DESC

   DESCRIPTION:
   This function fills the redundant components of a MATDATA_DESC.
   The functions 'ConstructMatOffsets' and 'SetScalMatSettings' are called.

   RETURN VALUE:
   INT
   .n    NUM_OK
   D*/
/****************************************************************************/

INT FillRedundantComponentsOfMD (MATDATA_DESC *md)
{
  ConstructMatOffsets(MD_ROWPTR(md),MD_COLPTR(md),MD_OFFSETPTR(md));
  SetCompactTypesOfMat(md);
  SetScalMatSettings(md);
  MD_SUCC_COMP(md) = MDCompsSubsequent(md);

  return (NUM_OK);
}

/****************************************************************************/
/*D
   GetFirstMatrix - return first matrix for a given multigrid

   SYNOPSIS:
   VECDATA_DESC *GetFirstMatrix (MULTIGRID *theMG)

   PARAMETERS:
   .  theMG - multigrid

   DESCRIPTION:
   This function returns the first matrix for a given multigrid

   RETURN VALUE:
   VECDATA_DESC *
   .n    NULL if there is none
   D*/
/****************************************************************************/

MATDATA_DESC *GetFirstMatrix (MULTIGRID *theMG)
{
  ENVITEM *item;

  if (ChangeEnvDir("/Multigrids") == NULL) return (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (NULL);
  item = (ENVITEM *)ChangeEnvDir("Matrices");
  if (item == NULL) return (NULL);
  for (item=ENVITEM_DOWN(item); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == MatrixVarID)
      return ((MATDATA_DESC *)item);

  return (NULL);
}

/****************************************************************************/
/*D
   GetNextMatrix - return next matrix for a given matrix

   SYNOPSIS:
   VECDATA_DESC *GetNextMatrix (VECDATA_DESC *vd)

   PARAMETERS:
   .  theMG - multigrid

   DESCRIPTION:
   This function returns the next matrix for a given matrix

   RETURN VALUE:
   VECDATA_DESC *
   .n    NULL if there is no successor
   D*/
/****************************************************************************/

MATDATA_DESC *GetNextMatrix (MATDATA_DESC *md)
{
  ENVITEM *item;

  for (item=NEXT_ENVITEM((ENVITEM *)md); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == MatrixVarID)
      return ((MATDATA_DESC *)item);

  return (NULL);
}

static INT GetNewMatrixName (MULTIGRID *theMG, char *name)
{
  MATDATA_DESC *md;
  char buffer[NAMESIZE];
  INT i;

  for (i=0; i<MAX_NAMES; i++) {
    sprintf(buffer,"mat%02d",i);
    for (md = GetFirstMatrix(theMG); md != NULL; md = GetNextMatrix(md))
      if (strcmp(ENVITEM_NAME(md),buffer) == 0) break;
    if (md == NULL) break;
  }
  if (i == MAX_NAMES) REP_ERR_RETURN (1);
  strcpy(name,buffer);

  return (NUM_OK);
}

/****************************************************************************/
/*D
   CreateMatDesc - create a matrix and fill extra data

   SYNOPSIS:
   MATDATA_DESC *CreateMatDesc (MULTIGRID *theMG, char *name, char *compNames,
   SHORT *RowsInType, SHORT *ColsInType);

   PARAMETERS:
   .  theMG - create matrix for this multigrid
   .  name - create matrix with this name
   .  RowsInType - 'MATDATA_DESC' specification
   .  ColsInType - 'MATDATA_DESC' specification
   .  compNames - (optional) matrix of component names (in the canonic type order)
               two char each (NULL pointer for no names)

   DESCRIPTION:
   This function creates a 'MATDATA_DESC' and fills its components.

   RETURN VALUE:
   MATDATA_DESC *
   .n    pointer to 'MATDATA_DESC' environment item
   .n     NULL if failed
   D*/
/****************************************************************************/

MATDATA_DESC *CreateMatDesc (MULTIGRID *theMG, const char *name, const char *compNames,
                             const SHORT *RowsInType, const SHORT *ColsInType)
{
  MATDATA_DESC *md;
  SHORT offset[NMATOFFSETS],*Comp;
  char buffer[NAMESIZE];
  INT i,j,tp,ncmp,size;

  if (theMG == NULL)
    REP_ERR_RETURN (NULL);

  if (ChangeEnvDir("/Multigrids") == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir("Matrices") == NULL) {
    MakeEnvItem("Matrices",MatrixDirID,sizeof(ENVDIR));
    if (ChangeEnvDir("Matrices") == NULL) REP_ERR_RETURN (NULL);
  }
  ConstructMatOffsets(RowsInType,ColsInType,offset);
  ncmp = offset[NMATTYPES];
  if (ncmp <= 0) REP_ERR_RETURN (NULL);
  size = sizeof(MATDATA_DESC)+(ncmp-1)*sizeof(SHORT);
  if (name != NULL)
    strcpy(buffer,name);
  else if (GetNewMatrixName(theMG,buffer)) REP_ERR_RETURN (NULL);
  md = (MATDATA_DESC *) MakeEnvItem (buffer,MatrixVarID,size);
  if (md == NULL) REP_ERR_RETURN (NULL);
  if (compNames==NULL)
    memcpy(VM_COMP_NAMEPTR(md),NoMatNames,2*MIN(ncmp,MAX_MAT_COMP));
  else
    memcpy(VM_COMP_NAMEPTR(md),compNames,2*MIN(ncmp,MAX_MAT_COMP));

  /* fill data in mat data desc */
  MD_MG(md) = theMG;
  i = 0;
  Comp = VM_COMPPTR(md);
  for (tp=0; tp<NMATTYPES; tp++) {
    ASSERT(RowsInType[tp]<A_REASONABLE_NUMBER);
    ASSERT(ColsInType[tp]<A_REASONABLE_NUMBER);
    MD_ROWS_IN_MTYPE(md,tp) = RowsInType[tp];
    MD_COLS_IN_MTYPE(md,tp) = ColsInType[tp];
    MD_MCMPPTR_OF_MTYPE(md,tp) = Comp + offset[tp];
    for (j=0; j<MAX_NDOF_MOD_32*32; j++) {
      if (i >= offset[tp+1]) break;
      if (j*sizeof(DOUBLE) >=
          FMT_S_MAT_TP(MGFORMAT(theMG),MatrixType[MTYPE_RT(tp)][MTYPE_CT(tp)]))
        REP_ERR_RETURN (NULL);
      if (READ_DR_MAT_FLAG(theMG,tp,j)) continue;
      Comp[i++] = j;
      SET_DR_MAT_FLAG(theMG,tp,j);
    }
  }
  for (tp=0; tp<NMATOFFSETS; tp++)
    MD_MTYPE_OFFSET(md,tp) = offset[tp];

  if (FillRedundantComponentsOfMD(md))
    return (NULL);
  VM_LOCKED(md) = 0;

  return (md);
}

/****************************************************************************/
/*D
   CreateSubMatDesc - create a matrix and fill extra data

   SYNOPSIS:
   MATDATA_DESC *CreateSubMatDesc (MULTIGRID *theMG, const char *name,
                                                                const SHORT *RowsInType, const SHORT *ColsInType,
                                                                const SHORT *Comps, const char *CompNames)

   PARAMETERS:
   .  theMG - create sub-matrix for this multigrid
   .  theMD - given matrix
   .  name - create sub-matrix with this name (maybe NULL for default name)
   .  RowsInType - 'MATDATA_DESC' specification
   .  ColsInType - 'MATDATA_DESC' specification
   .  Comps - the Comps of theMD which make up the sub-matdesc
   .  compNames - (optional) matrix of component names (in the canonic type order)
               two char each (NULL pointer for no names)

   DESCRIPTION:
   This function creates a 'MATDATA_DESC' and fills its components.

   RETURN VALUE:
   MATDATA_DESC *
   .n    pointer to 'MATDATA_DESC' environment item
   .n     NULL if failed
   D*/
/****************************************************************************/

MATDATA_DESC *CreateSubMatDesc (MULTIGRID *theMG, const char *name,
                                const SHORT *RowsInType, const SHORT *ColsInType,
                                const SHORT *Comps, const char *CompNames)
{
  MATDATA_DESC *md;
  SHORT offset[NMATOFFSETS];
  INT k,j,tp,ncmp,size;
  char buffer[NAMESIZE];

  if (theMG == NULL)
    REP_ERR_RETURN (NULL);

  if (ChangeEnvDir("/Multigrids") == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir("Matrices") == NULL) REP_ERR_RETURN (NULL);
  ConstructMatOffsets(RowsInType,ColsInType,offset);
  ncmp = offset[NMATTYPES];
  if (ncmp <= 0) REP_ERR_RETURN (NULL);
  size = sizeof(MATDATA_DESC)+(ncmp-1)*sizeof(SHORT);
  if (name != NULL)
    strcpy(buffer,name);
  else if (GetNewMatrixName(theMG,buffer)) REP_ERR_RETURN (NULL);
  md = (MATDATA_DESC *) MakeEnvItem (buffer,MatrixVarID,size);
  if (md == NULL) REP_ERR_RETURN (NULL);

  /* fill data in mat data desc */
  MD_MG(md) = theMG;
  strncpy(VM_COMP_NAMEPTR(md),CompNames,2*ncmp);
  k = 0;
  for (tp=0; tp<NMATTYPES; tp++) {
    MD_ROWS_IN_MTYPE(md,tp) = RowsInType[tp];
    MD_COLS_IN_MTYPE(md,tp) = ColsInType[tp];
    MD_MCMPPTR_OF_MTYPE(md,tp) = VM_COMPPTR(md) + offset[tp];
    for (j=0; j<RowsInType[tp]*ColsInType[tp]; j++) {
      MD_MCMP_OF_MTYPE(md,tp,j) = Comps[k++];
    }
  }
  for (tp=0; tp<NMATOFFSETS; tp++)
    MD_MTYPE_OFFSET(md,tp) = offset[tp];

  if (FillRedundantComponentsOfMD(md))
    return (NULL);
  VM_LOCKED(md) = 0;

  return (md);
}

/****************************************************************************/
/*D
   DisplayMatDataDesc - Display MATDATA_DESC entries

   SYNOPSIS:
   INT DisplayMatDataDesc (const MATDATA_DESC *md, char *buffer)

   PARAMETERS:
   .  fmt - associated format for names of abstract types
   .  md - MATDATA_DESC to display
   .  buffer - print here

   DESCRIPTION:
   This function displays the entries of a MATDATA_DESC: comp-names, comp-positions etc.

   RETURN VALUE:
   INT
   .n      0: ok
   .n      else: error
   D*/
/****************************************************************************/

INT DisplayMatDataDesc (const MATDATA_DESC *md, char *buffer)
{
  const FORMAT *fmt;
  const SHORT *offset;
  const char *cn;
  INT rt,ct,mtp,i,j,nc,maxr[NVECTYPES],maxc[NVECTYPES];

  if (md==NULL) REP_ERR_RETURN (1);

  buffer += sprintf(buffer,"contents of matrix symbol '%s'\n",ENVITEM_NAME(md));

  fmt = MGFORMAT(MD_MG(md));
  cn = VM_COMP_NAMEPTR(md);
  offset = MD_OFFSETPTR(md);
  if (cn[0]==' ')
    cn = NULL;
  else
    for (i=0; i<offset[NMATTYPES]; i++)
      if (cn[i]=='\0') {
        cn = NULL;
        break;
      }

  for (rt=0; rt<NVECTYPES; rt++)
  {
    maxr[rt] = 0;
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
        maxr[rt] = MAX(maxr[rt],MD_ROWS_IN_RT_CT(md,rt,ct));
  }

  /* headline for col types */
  buffer += sprintf(buffer,"  ");
  for (ct=0; ct<NVECTYPES; ct++)
  {
    maxc[ct] = 0;
    for (rt=0; rt<NVECTYPES; rt++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
        maxc[ct] = MAX(maxc[ct],MD_COLS_IN_RT_CT(md,rt,ct));
    for (j=0; j<maxc[ct]; j++)
      buffer += sprintf(buffer," %s%c ",(j) ? "" : "|",(j) ? ' ' : FMT_VTYPE_NAME(fmt,ct));
  }
  buffer += sprintf(buffer,"\n--");
  for (ct=0; ct<NVECTYPES; ct++)
    for (j=0; j<maxc[ct]; j++)
      buffer += sprintf(buffer,"-%s--",(j) ? "" : "-");

  for (rt=0; rt<NVECTYPES; rt++)
  {
    for (i=0; i<maxr[rt]; i++)
    {
      /* compname line */
      buffer += sprintf(buffer,"\n%c ",(i) ? ' ' : FMT_VTYPE_NAME(fmt,rt));
      if (cn!=NULL)
      {
        for (ct=0; ct<NVECTYPES; ct++)
        {
          j = 0;
          if (MD_ISDEF_IN_RT_CT(md,rt,ct))
          {
            mtp = MTP(rt,ct);
            nc  = MD_COLS_IN_MTYPE(md,mtp);
            for (; j<nc; j++)
              buffer += sprintf(buffer," %s%c%c",(j) ? "" : "|",cn[2*(offset[mtp]+i*nc+j)],cn[2*(offset[mtp]+i*nc+j)+1]);
          }
          for (; j<maxc[ct]; j++)
            buffer += sprintf(buffer," %s  ",(j) ? "" : "|");
        }
        buffer += sprintf(buffer,"\n  ");
      }
      /* comp position line */
      for (ct=0; ct<NVECTYPES; ct++)
      {
        j = 0;
        if (MD_ISDEF_IN_RT_CT(md,rt,ct))
        {
          mtp = MTP(rt,ct);
          for (; j<MD_COLS_IN_MTYPE(md,mtp); j++)
            buffer += sprintf(buffer," %s%2d",(j) ? "" : "|",MD_IJ_CMP_OF_MTYPE(md,mtp,i,j));
        }
        for (; j<maxc[ct]; j++)
          buffer += sprintf(buffer," %s  ",(j) ? "" : "|");
      }
    }
    if (maxr[rt]>0)
    {
      /* type seperator line */
      buffer += sprintf(buffer,"\n--");
      for (ct=0; ct<NVECTYPES; ct++)
        for (j=0; j<maxc[ct]; j++)
          buffer += sprintf(buffer,"-%s--",(j) ? "" : "-");
    }
  }
  buffer += sprintf(buffer,"\n");

  if (MD_IS_SCALAR(md))
  {
    buffer += sprintf(buffer,"\nmatsym is scalar:\n");
    buffer += sprintf(buffer,"  comp %2d\n",MD_SCALCMP(md));
    buffer += sprintf(buffer,"  rmsk %2d\n",MD_SCAL_RTYPEMASK(md));
    buffer += sprintf(buffer,"  cmsk %2d\n",MD_SCAL_CTYPEMASK(md));
  }

  buffer += sprintf(buffer,"\n");

  return (NUM_OK);
}

/****************************************************************************/
/*D
   GetMatDataDescByName - find matrix data descriptor

   SYNOPSIS:
   MATDATA_DESC *GetMatDataDescByName (MULTIGRID *theMG, char *name);

   PARAMETERS:
   .  theMG - create vector for this multigrid
   .  name - name of a matrix

   DESCRIPTION:
   This function finds a matrix by name.

   RETURN VALUE:
   MATDATA_DESC *
   .n      pointer to the matrix
   .n      NULL if there is no matrix of this name in the multigrid
   D*/
/****************************************************************************/

MATDATA_DESC *GetMatDataDescByName (const MULTIGRID *theMG, char *name)
{
  if (ChangeEnvDir("/Multigrids") == NULL) return (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (NULL);
  return((MATDATA_DESC *) SearchEnv(name,"Matrices",
                                    MatrixVarID,MatrixDirID));
}

/****************************************************************************/
/*D
   AllocMDFromVD - dynamic matrix allocation

   SYNOPSIS:
   INT AllocMDFromVD (MULTIGRID *theMG, INT fl, INT tl,
   VECDATA_DESC *x, VECDATA_DESC *y, MATDATA_DESC **new_desc);

   PARAMETERS:
   .  theMG - create vector for this multigrid
   .  fl - from level
   .  tl - to level
   .  x - template vector for row components
   .  y - template vector for column components
   .  new_desc - new matrix

   DESCRIPTION:
   This function allocates a new matrix.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

static INT CompMatDesc (const MATDATA_DESC *md,
                        const SHORT *RowsInType, const SHORT *ColsInType)
{
  INT tp;

  for (tp=0; tp<NMATTYPES; tp++) {
    if (MD_COLS_IN_MTYPE(md,tp) != ColsInType[tp])
      return(1);
    if (MD_ROWS_IN_MTYPE(md,tp) != RowsInType[tp])
      return(1);
  }

  return(0);
}

static INT AllocMatDesc (MULTIGRID *theMG, INT fl, INT tl, MATDATA_DESC *md)
{
  GRID *theGrid;
  INT i,j,tp;

  if (md == NULL) return(1);
  for (i=fl; i<=tl; i++) {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (tp=0; tp<NMATTYPES; tp++)
      for (j=0; j<MD_NCMPS_IN_MTYPE(md,tp); j++)
        if (READ_DR_MAT_FLAG(theGrid,tp,MD_MCMP_OF_MTYPE(md,tp,j)))
          return(1);
  }
  for (i=fl; i<=tl; i++) {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (tp=0; tp<NMATTYPES; tp++)
      for (j=0; j<MD_NCMPS_IN_MTYPE(md,tp); j++)
        SET_DR_MAT_FLAG(theGrid,tp,MD_MCMP_OF_MTYPE(md,tp,j));
  }

  return(0);
}

/****************************************************************************/
/*D
   AllocMDFromMRowMCol - dynamically allocate matrix descriptor from given RowsInType, ColsInType and compNames
                                in mtype style

   SYNOPSIS:
   INT AllocMDFromMRowMCol (MULTIGRID *theMG, INT fl, INT tl,
                                   const SHORT *RowsInType,const SHORT *ColsInType,const char *compNames,MATDATA_DESC **new_desc)

   PARAMETERS:
   .  theMG - allocate descriptor for this multigrid
   .  fl - from level
   .  tl - to level
   .  RowsInType - rows per mtype
   .  ColsInType - cols per mtype
   .  compNames  - names of comps (may be NULL)
   .  new_desc   - handle for descriptor

   DESCRIPTION:
   This function allocates a matrix descriptor from
   given RowsInType, ColsInType and compNames. If there is no existing
   decriptor a new one is created.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT AllocMDFromMRowMCol (MULTIGRID *theMG, INT fl, INT tl,
                         const SHORT *RowsInType,const SHORT *ColsInType,const char *compNames,MATDATA_DESC **new_desc)
{
  MATDATA_DESC *md;

  /* nothing to do if it is locked */
  if (*new_desc != NULL)
    if (VM_LOCKED(*new_desc))
      return (NUM_OK);

  /* is there a freed Matrix we can use */
  if (AllocMatDesc(theMG,fl,tl,*new_desc)) {
    for (md = GetFirstMatrix(theMG); md != NULL; md = GetNextMatrix(md)) {
      if (VM_LOCKED(md)) continue;
      if (CompMatDesc(md,RowsInType,
                      ColsInType))
        continue;
      if (!AllocMatDesc(theMG,fl,tl,md)) {
        *new_desc = md;
        return (NUM_OK);
      }
    }
    /* create a new Matrix */
    *new_desc = CreateMatDesc(theMG,NULL,compNames,
                              RowsInType,
                              ColsInType);
    if (*new_desc == NULL) {
      PrintErrorMessage('E',"AllocMDFromMRowMCol","cannot create MatDesc\n");
      REP_ERR_RETURN(1);
    }
    if (AllocMatDesc(theMG,fl,tl,*new_desc)) {
      PrintErrorMessage('E',"AllocMDFromMRowMCol","cannot allocate MatDesc\n");
      REP_ERR_RETURN(1);
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   AllocMDFromVRowVCol - dynamically allocate matrix descriptor from given RowsInType, ColsInType in vtype style

   SYNOPSIS:
   INT AllocMDFromVRowVCol (MULTIGRID *theMG, INT fl, INT tl,
                                                   const SHORT *RowsInType,const SHORT *ColsInType,MATDATA_DESC **new_desc)

   PARAMETERS:
   .  theMG - allocate descriptor for this multigrid
   .  fl - from level
   .  tl - to level
   .  RowsInType - rows per vtype
   .  ColsInType - cols per vtype
   .  compNames  - names of comps (may be NULL)
   .  new_desc   - handle for descriptor

   DESCRIPTION:
   This function allocates a matrix descriptor from
   given RowsInType, ColsInType as if they where from vectors the matrix is operating
   from and to. If there is no existing decriptor a new one is created.

   SEE ALSO:
   AllocMDFromMRowMCol

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT AllocMDFromVRowVCol (MULTIGRID *theMG, INT fl, INT tl,
                         const SHORT *RowsInType,const SHORT *ColsInType,MATDATA_DESC **new_desc)
{
  return (AllocMDFromMRowMCol(theMG,fl,tl,RowsInType,ColsInType,NULL,new_desc));
}

/****************************************************************************/
/*D
   AllocMDFromVD - dynamically allocate matrix descriptor from given vector descriptors

   SYNOPSIS:
   INT AllocMDFromVD (MULTIGRID *theMG, INT fl, INT tl,
                                   const VECDATA_DESC *x, const VECDATA_DESC *y, MATDATA_DESC **new_desc)

   PARAMETERS:
   .  theMG - allocate descriptor for this multigrid
   .  fl - from level
   .  tl - to level
   .  x - source vector the matrix is operating on
   .  y - destination vector the matrix is operating to
   .  new_desc   - handle for descriptor

   DESCRIPTION:
   This function allocates a matrix descriptor from
   given vector descriptors the matrix is operating
   from and to. If there is no existing decriptor a new one is created.

   SEE ALSO:
   AllocMDFromMRowMCol

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT AllocMDFromVD (MULTIGRID *theMG, INT fl, INT tl,
                   const VECDATA_DESC *x, const VECDATA_DESC *y, MATDATA_DESC **new_desc)
{
  INT i,j,tp;
  SHORT RowsInType[NMATTYPES];
  SHORT ColsInType[NMATTYPES];

  /* nothing to do if it is locked */
  if (*new_desc != NULL)
    if (VM_LOCKED(*new_desc))
      return (NUM_OK);

  /* compute RowsInType, ColsInType */
  if (AllocMatDesc(theMG,fl,tl,*new_desc)==NUM_OK)
    return (NUM_OK);

  /* translate vd-x and vd-y into RowsInType and ColsInType */
  for (i=0; i<NVECTYPES; i++)
    for (j=0; j<NVECTYPES; j++) {
      tp = MTP(i,j);
      if (VD_NCMPS_IN_TYPE(x,i)*VD_NCMPS_IN_TYPE(y,j) > 0) {
        RowsInType[tp] = VD_NCMPS_IN_TYPE(x,i);
        ColsInType[tp] = VD_NCMPS_IN_TYPE(y,j);
      }
      else
        RowsInType[tp] = ColsInType[tp] = 0;
    }

  return (AllocMDFromMRowMCol(theMG,fl,tl,RowsInType,ColsInType,NULL,new_desc));
}

/****************************************************************************/
/*D
   AllocMDFromMD - dynamically allocate matrix descriptor from given matrix descriptor

   SYNOPSIS:
   INT AllocMDFromMD (MULTIGRID *theMG, INT fl, INT tl,
                                   const MATDATA_DESC *md, MATDATA_DESC **new_desc)

   PARAMETERS:
   .  theMG - allocate descriptor for this multigrid
   .  fl - from level
   .  tl - to level
   .  md - tempplate matrix descriptor
   .  new_desc   - handle for descriptor

   DESCRIPTION:
   This function allocates a matrix descriptor from
   a given matrix descriptor. If there is no existing decriptor a new one is created.

   SEE ALSO:
   AllocMDFromMRowMCol

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT AllocMDFromMD (MULTIGRID *theMG, INT fl, INT tl,
                   const MATDATA_DESC *md, MATDATA_DESC **new_desc)
{
  return (AllocMDFromMRowMCol(theMG,fl,tl,
                              md->RowsInType,
                              md->ColsInType,
                              md->compNames,
                              new_desc));
}

/****************************************************************************/
/*D
   LockMD - protect matrix against removal or deallocation

   SYNOPSIS:
   INT LockMD (MATDATA_DESC *md)

   PARAMETERS:
   .  md - matrix descriptor

   DESCRIPTION:
   This function protects a matrix against removal or deallocation.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT LockMD (MATDATA_DESC *md)
{
  VM_LOCKED(md) = VM_IS_LOCKED;
  return (0);
}

/****************************************************************************/
/*D
   TransmitLockStatusMD - ...

   SYNOPSIS:
   INT TransmitLockStatusMD (const MATDATA_DESC *md, MATDATA_DESC *smd)

   PARAMETERS:
   .  md  - matrix descriptor
   .  smd - sub matrix descriptor

   DESCRIPTION:
   This function ...

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT TransmitLockStatusMD (const MATDATA_DESC *md, MATDATA_DESC *smd)
{
  if (!VM_LOCKED(md) && VM_LOCKED(smd))
    REP_ERR_RETURN(1);
  VM_LOCKED(smd) = VM_LOCKED(md);

  return (0);
}

/****************************************************************************/
/*D
   FreeMD - dynamic matrix deallocation

   SYNOPSIS:
   INT FreeMD (MULTIGRID *theMG, INT fl, INT tl, MATDATA_DESC *md);

   PARAMETERS:
   .  theMG - free vector for this multigrid
   .  fl - from level
   .  tl - to level
   .  md - matrix descriptor

   DESCRIPTION:
   This function deallocates a matrix.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT FreeMD (MULTIGRID *theMG, INT fl, INT tl, MATDATA_DESC *md)
{
  GRID *theGrid;
  INT i,j,tp;

  if (md==NULL) return (NUM_OK);
  if (VM_LOCKED(md)) return (NUM_OK);
  for (i=fl; i<=tl; i++) {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (tp=0; tp<NMATTYPES; tp++)
      for (j=0; j<MD_NCMPS_IN_MTYPE(md,tp); j++)
        CLEAR_DR_MAT_FLAG(theGrid,tp,MD_MCMP_OF_MTYPE(md,tp,j));
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   DisposeMD - remove MATDATA_DESC from objects of multigrid

   SYNOPSIS:
   INT DisposeMD (MATDATA_DESC *md)

   PARAMETERS:
   .  md - matrix descriptor

   DESCRIPTION:
   This function removes a matrix descriptor from the objects of the multigrid. The part
   of environment memory it occupies is freed. The functions calls
   'RemoveEnvItem'.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occurred
 */
/****************************************************************************/

INT DisposeMD (MATDATA_DESC *md)
{
  if (md==NULL) REP_ERR_RETURN (NUM_ERROR);
  if (VM_LOCKED(md)) REP_ERR_RETURN (NUM_ERROR);

  ENVITEM_LOCKED(md) = 0;

  if (ChangeEnvDir("/Multigrids") == NULL) REP_ERR_RETURN (0);
  if (ChangeEnvDir(ENVITEM_NAME(VD_MG(md))) == NULL) REP_ERR_RETURN (0);
  if (ChangeEnvDir("Matrices") == NULL) REP_ERR_RETURN (0);

  if (RemoveEnvItem((ENVITEM*)md))
    REP_ERR_RETURN (0);

  return (0);
}

/****************************************************************************/
/*D
        MDinterfaceDesc - an interface MATDATA_DESC is created

        SYNOPSIS:
        INT MDinterfaceDesc (const MATDATA_DESC *md, const MATDATA_DESC *mds, MATDATA_DESC **mdi)

    PARAMETERS:
   .   md			- make a sub desc of this MATDATA_DESC
   .   mds			- an existing sub desc of md
   .   mdi			- handle to new interface desc

        DESCRIPTION:
        This function creates a sub descriptor to a given MATDATA_DESC such that all components
        from mds are taken of types in which mds is defined but where mds is a true subset of md.

        RETURN VALUE:
        INT
   .n   0: ok
   .n      n: if an error occured
   D*/
/****************************************************************************/

INT MDinterfaceDesc (const MATDATA_DESC *md, const MATDATA_DESC *mds, MATDATA_DESC **mdi)
{
  SHORT SubComp[MAX_MAT_COMP],SubRCmp[NMATTYPES],SubCCmp[NMATTYPES];
  INT i,k,l,n,ns,tp;
  char SubName[2*MAX_MAT_COMP],buffer[NAMESIZE];

  /* generate name and see if desc already exists */
  strcpy(buffer,ENVITEM_NAME(mds));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,"i");
  *mdi = GetMatDataDescByName(MD_MG(md),buffer);
  if (*mdi != NULL) {
    if (TransmitLockStatusMD(mds,*mdi))
      REP_ERR_RETURN(1);
    return(0);
  }

  k = 0;
  for (tp=0; tp<NMATTYPES; tp++)
  {
    SubRCmp[tp] = SubCCmp[tp] = 0;
    if (MD_ISDEF_IN_MTYPE(mds,tp))
    {
      if (!MD_ISDEF_IN_MTYPE(md,tp))
        REP_ERR_RETURN (1);

      n  = MD_NCMPS_IN_MTYPE(md,tp);
      ns = MD_NCMPS_IN_MTYPE(mds,tp);
      if (ns<n)
      {
        /* copy all components from mds */
        for (i=0; i<ns; i++)
        {
          ASSERT(k<MAX_MAT_COMP);

          SubComp[k] = MD_MCMP_OF_MTYPE(mds,tp,i);
          l = MD_MTYPE_OFFSET(mds,tp)+i;
          SubName[2*k]   = VM_COMP_NAME(mds,2*l);
          SubName[2*k+1] = VM_COMP_NAME(mds,2*l+1);
          k++;
        }
        SubRCmp[tp] = MD_ROWS_IN_MTYPE(mds,tp);
        SubCCmp[tp] = MD_COLS_IN_MTYPE(mds,tp);
      }
      else if (ns!=n)
        /* md does not contain mds */
        REP_ERR_RETURN (1);
    }
  }

  *mdi = CreateSubMatDesc(MD_MG(md),buffer,SubRCmp,SubCCmp,SubComp,SubName);
  if (*mdi == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusMD(md,*mdi))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/*D
        MDinterfaceCoCoupleDesc - an interface MATDATA_DESC is created describing couplings

        SYNOPSIS:
        INT MDinterfaceCoCoupleDesc (const MATDATA_DESC *md, const MATDATA_DESC *mds, MATDATA_DESC **mdi)

    PARAMETERS:
   .   md			- make a sub desc of this MATDATA_DESC
   .   mds			- an existing sub desc of md
   .   mdi			- handle to new coupling desc

        DESCRIPTION:
        This function creates a sub descriptor to a given MATDATA_DESC such that all components
        from md which are not in mds but in the same rows as mds will be described.
        Components are taken only in types in which mds is defined but where mds is a true subset of md.

        RETURN VALUE:
        INT
   .n   0: ok
   .n      n: if an error occured
   D*/
/****************************************************************************/

INT MDinterfaceCoCoupleDesc (const MATDATA_DESC *md, const MATDATA_DESC *mds, MATDATA_DESC **mdi)
{
  SHORT SubComp[MAX_MAT_COMP],SubRCmp[NMATTYPES],SubCCmp[NMATTYPES];
  INT i,j,jj,k,l,n,ns,nr,co_nr,nc,rt,ct,cmp;
  INT RowUsed[MAX_VEC_COMP];
  char SubName[2*MAX_MAT_COMP],buffer[NAMESIZE];

  /* generate name and see if desc already exists */
  strcpy(buffer,ENVITEM_NAME(mds));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,"icc");
  *mdi = GetMatDataDescByName(MD_MG(md),buffer);
  if (*mdi != NULL) {
    if (TransmitLockStatusMD(mds,*mdi))
      REP_ERR_RETURN(1);
    return(0);
  }

  k = 0;
  for (rt=0; rt<NVECTYPES; rt++)
  {
    for (ct=0; ct<NVECTYPES; ct++)
      SubRCmp[MTP(rt,ct)] = SubCCmp[MTP(rt,ct)] = 0;

    /* find type rows used my mds */
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(mds,rt,ct))
        break;
    if (ct>=NVECTYPES)
      continue;

    /* find rows in (rt,ct) used by mds relative to md */
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(mds,rt,ct))
      {
        if (!MD_ISDEF_IN_RT_CT(md,rt,ct))
          REP_ERR_RETURN(1);

        n  = MD_NCMPS_IN_RT_CT(md,rt,ct);
        ns = MD_NCMPS_IN_RT_CT(mds,rt,ct);
        if (n==ns)
          continue;
        if (ns>n)
          REP_ERR_RETURN(1);

        nr = MD_ROWS_IN_RT_CT(md,rt,ct);
        nc = MD_COLS_IN_RT_CT(md,rt,ct);
        for (i=0; i<nr; i++)
          RowUsed[i] = FALSE;
        for (jj=0; jj<ns; jj++)
        {
          cmp = MD_MCMP_OF_RT_CT(mds,rt,ct,jj);
          for (i=0; i<nr; i++)
            for (j=0; j<nc; j++)
              if (cmp==MD_IJ_CMP_OF_RT_CT(md,rt,ct,i,j))
                RowUsed[i] = TRUE;
        }
      }
    co_nr = 0;
    for (i=0; i<nr; i++)
      if (RowUsed[i])
        co_nr++;

    /* now copy all md-comps from used rows not in mds */
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
      {
        nr = MD_ROWS_IN_RT_CT(md,rt,ct);
        nc = MD_COLS_IN_RT_CT(md,rt,ct);
        if (!MD_ISDEF_IN_RT_CT(mds,rt,ct))
        {
          /* copy all comps of used rows */
          for (i=0; i<nr; i++)
            if (RowUsed[i])
              for (j=0; j<nc; j++)
              {
                ASSERT(k<MAX_MAT_COMP);

                SubComp[k] = MD_IJ_CMP_OF_RT_CT(md,rt,ct,i,j);
                l = MD_MTYPE_OFFSET(md,MTP(rt,ct))+i*nr+j;
                SubName[2*k]   = VM_COMP_NAME(md,2*l);
                SubName[2*k+1] = VM_COMP_NAME(md,2*l+1);
                k++;
              }
          SubRCmp[MTP(rt,ct)] = co_nr;
          SubCCmp[MTP(rt,ct)] = nc;
        }
        else
        {
          if (nc==MD_COLS_IN_RT_CT(mds,rt,ct))
            /* no comps to copy */
            continue;

          /* copy all comps of used rows not in mds */
          ns = MD_NCMPS_IN_RT_CT(mds,rt,ct);
          for (i=0; i<nr; i++)
            if (RowUsed[i])
              for (j=0; j<nc; j++)
              {
                cmp = MD_IJ_CMP_OF_RT_CT(md,rt,ct,i,j);

                for (jj=0; jj<ns; jj++)
                  if (cmp==MD_MCMP_OF_RT_CT(mds,rt,ct,jj))
                    break;
                if (jj<ns)
                  continue;

                ASSERT(k<MAX_MAT_COMP);

                SubComp[k] = cmp;
                l = MD_MTYPE_OFFSET(md,MTP(rt,ct))+i*nr+j;
                SubName[2*k]   = VM_COMP_NAME(md,2*l);
                SubName[2*k+1] = VM_COMP_NAME(md,2*l+1);
                k++;
              }
          SubRCmp[MTP(rt,ct)] = co_nr;
          SubCCmp[MTP(rt,ct)] = nc-MD_COLS_IN_RT_CT(mds,rt,ct);
        }
      }
  }

  *mdi = CreateSubMatDesc(MD_MG(md),buffer,SubRCmp,SubCCmp,SubComp,SubName);
  if (*mdi == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusMD(md,*mdi))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/*D
   MD_rows_in_ro_co - return number of row comps in row/col object if unique, -1 else

   SYNOPSIS:
   INT MD_rows_in_ro_co (const MATDATA_DESC *md, INT rowobj, INT colobj, INT mode)

   PARAMETERS:
   .  md - data decsriptor
   .  rowobj - row object type
   .  colobj - col object type
   .  mode - STRICT or NON_STRICT

   DESCRIPTION:
   This function checks whether the number of row components described in 'rowobj'/'colobj'
   is the same for all mtypes using objects of 'rowobj'/'colobj' and returns it.
   If the number is not unique a -1 is returned.
   If not the whole domain is covered, a -2 is returned.
   The uniqueness of comps is not checked here.
   If mode is STRICT vectors are required in ALL vobjects of otype.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'otype'
   are not defined at all!

   RETURN VALUE:
   INT
   .n      number of components in objects of 'otype'
   .n      -1 if not unique
   .n      -2 if not the whole domain is covered
   D*/
/****************************************************************************/

INT MD_rows_in_ro_co_mod (const MATDATA_DESC *md, INT rowobj, INT colobj, INT mode)
{
  FORMAT *fmt;
  INT rt,ct,rot,cot,nrow,src_parts,dst_parts,i,n;

  fmt = MGFORMAT(MD_MG(md));
  rot = 1<<rowobj;
  cot = 1<<colobj;
  nrow = src_parts = dst_parts = 0;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
        if ((rot & FMT_T2O(fmt,rt)) && (cot & FMT_T2O(fmt,ct)))
        {
          if (nrow==0)
            nrow = MD_ROWS_IN_RT_CT(md,rt,ct);
          else
          if (MD_ROWS_IN_RT_CT(md,rt,ct)!=nrow)
            REP_ERR_RETURN (-1);
          src_parts |= FMT_T2P(fmt,rt);
          dst_parts |= FMT_T2P(fmt,ct);
        }

  if (mode==STRICT)
  {
    /* now check whether all parts are covered */
    src_parts &= dst_parts;
    n = BVPD_NPARTS(MG_BVPD(VD_MG(md)));
    for (i=0; i<n; i++)
      if (!(src_parts & (1<<i)))
        REP_ERR_RETURN (-2);
  }
  else if (mode!=NON_STRICT)
    REP_ERR_RETURN (1);

  return (nrow);
}

/****************************************************************************/
/*D
   MD_cols_in_ro_co - return number of col comps in row/col object if unique, -1 else

   SYNOPSIS:
   INT MD_cols_in_ro_co (const MATDATA_DESC *md, INT rowobj, INT colobj, INT mode)

   PARAMETERS:
   .  md - data decsriptor
   .  rowobj - row object type
   .  colobj - col object type
   .  mode - STRICT or NON_STRICT

   DESCRIPTION:
   This function checks whether the number of col components described in 'rowobj'/'colobj'
   is the same for all mtypes using objects of 'rowobj'/'colobj' and returns it.
   If the number is not unique a -1 is returned.
   If not the whole domain is covered, a -2 is returned.
   The uniqueness of comps is not checked here.
   If mode is STRICT vectors are required in ALL vobjects of otype.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'otype'
   are not defined at all!

   RETURN VALUE:
   INT
   .n      number of components in objects of 'otype'
   .n      -1 if not unique
   .n      -2 if not the whole domain is covered
   D*/
/****************************************************************************/

INT MD_cols_in_ro_co_mod (const MATDATA_DESC *md, INT rowobj, INT colobj, INT mode)
{
  FORMAT *fmt;
  INT rt,ct,rot,cot,ncol,src_parts,dst_parts,i,n;

  fmt = MGFORMAT(MD_MG(md));
  rot = 1<<rowobj;
  cot = 1<<colobj;
  ncol = src_parts = dst_parts = 0;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
        if ((rot & FMT_T2O(fmt,rt)) && (cot & FMT_T2O(fmt,ct)))
        {
          if (ncol==0)
            ncol = MD_COLS_IN_RT_CT(md,rt,ct);
          else
          if (MD_COLS_IN_RT_CT(md,rt,ct)!=ncol)
            REP_ERR_RETURN (-1);
          src_parts |= FMT_T2P(fmt,rt);
          dst_parts |= FMT_T2P(fmt,ct);
        }

  if (mode==STRICT)
  {
    /* now check whether all parts are covered */
    src_parts &= dst_parts;
    n = BVPD_NPARTS(MG_BVPD(VD_MG(md)));
    for (i=0; i<n; i++)
      if (!(src_parts & (1<<i)))
        REP_ERR_RETURN (-2);
  }
  else if (mode!=NON_STRICT)
    REP_ERR_RETURN (1);

  return (ncol);
}

/****************************************************************************/
/*D
   MD_rows_cols_in_ro_co - return number of row and col comps in row/col object if unique, -1 else

   SYNOPSIS:
   INT MD_rows_cols_in_ro_co (const MATDATA_DESC *md, INT rowobj, INT colobj, INT *nr, INT *nc, INT mode)

   PARAMETERS:
   .  md - data decsriptor
   .  rowobj - row object type
   .  colobj - col object type
   .  nr - number of rows if unique, unchanged if error
   .  nc - number of cols if unique, unchanged if error
   .  mode - STRICT or NON_STRICT

   DESCRIPTION:
   This function checks whether the number of row and col components described in 'rowobj'/'colobj'
   is the same for all mtypes using objects of 'rowobj'/'colobj'.
   If the number is not unique 1 is returned.
   If not the whole domain is covered, 2 is returned.
   The uniqueness of comps is not checked here.
   If mode is STRICT vectors are required in ALL vobjects of otype.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'otype'
   are not defined at all!

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if not unique
   .n      2 if not the whole domain is covered
   D*/
/****************************************************************************/

INT MD_rows_cols_in_ro_co_mod (const MATDATA_DESC *md, INT rowobj, INT colobj, INT *nr, INT *nc, INT mode)
{
  FORMAT *fmt;
  INT rt,ct,rot,cot,nrow,ncol,src_parts,dst_parts,i,n;

  fmt = MGFORMAT(MD_MG(md));
  rot = 1<<rowobj;
  cot = 1<<colobj;
  nrow = ncol = src_parts = dst_parts = 0;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
        if ((rot & FMT_T2O(fmt,rt)) && (cot & FMT_T2O(fmt,ct)))
        {
          if (nrow==0)
          {
            nrow = MD_ROWS_IN_RT_CT(md,rt,ct);
            ncol = MD_COLS_IN_RT_CT(md,rt,ct);
          }
          else
          {
            if (MD_ROWS_IN_RT_CT(md,rt,ct)!=nrow)
              REP_ERR_RETURN (1);
            if (MD_COLS_IN_RT_CT(md,rt,ct)!=ncol)
              REP_ERR_RETURN (1);
          }
          src_parts |= FMT_T2P(fmt,rt);
          dst_parts |= FMT_T2P(fmt,ct);
        }

  if (mode==STRICT)
  {
    /* now check whether all parts are covered */
    src_parts &= dst_parts;
    n = BVPD_NPARTS(MG_BVPD(VD_MG(md)));
    for (i=0; i<n; i++)
      if (!(src_parts & (1<<i)))
        REP_ERR_RETURN (2);
  }
  else if (mode!=NON_STRICT)
    REP_ERR_RETURN (1);

  *nr = nrow;
  *nc = ncol;

  return (0);
}

/****************************************************************************/
/*D
   MD_mcmp_of_ro_co - return comp in row/col object if unique, -1 else

   SYNOPSIS:
   INT MD_mcmp_of_ro_co (const MATDATA_DESC *md, INT rowobj, INT colobj, INT i, INT mode)

   PARAMETERS:
   .  md - data decsriptor
   .  rowobj - row object type
   .  colobj - col object type
   .  i - component number
   .  mode - STRICT or NON_STRICT

   DESCRIPTION:
   This function checks whether the offset of component number i described in 'rowobj'/'colobj'
   is the same for all mtypes using objects of 'rowobj'/'colobj' and returns it.
   If the offset is not unique a -1 is returned.
   If not the whole domain is covered, a -2 is returned.
   If mode is STRICT vectors are required in ALL vobjects of otype.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'otype'
   are not defined at all!

   RETURN VALUE:
   INT
   .n      number of components in objects of 'otype'
   .n      -1 if not unique
   .n      -2 if not the whole domain is covered
   D*/
/****************************************************************************/

INT MD_mcmp_of_ro_co_mod (const MATDATA_DESC *md, INT rowobj, INT colobj, INT i, INT mode)
{
  FORMAT *fmt;
  INT rt,ct,off,rot,cot,nrow,ncol,src_parts,dst_parts,j,n;

  fmt = MGFORMAT(MD_MG(md));
  rot = 1<<rowobj;
  cot = 1<<colobj;
  off = nrow = ncol = src_parts = dst_parts = 0;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
        if ((rot & FMT_T2O(fmt,rt)) && (cot & FMT_T2O(fmt,ct)))
        {
          if (nrow==0)
          {
            nrow = MD_ROWS_IN_RT_CT(md,rt,ct);
            ncol = MD_COLS_IN_RT_CT(md,rt,ct);
            off  = MD_MCMP_OF_RT_CT(md,rt,ct,i);
            if (i>=nrow*ncol)
              REP_ERR_RETURN (-1);
          }
          else
          {
            if (MD_ROWS_IN_RT_CT(md,rt,ct)!=nrow)
              REP_ERR_RETURN (-1);
            if (MD_COLS_IN_RT_CT(md,rt,ct)!=ncol)
              REP_ERR_RETURN (-1);
            if (MD_MCMP_OF_RT_CT(md,rt,ct,i)!=off)
              REP_ERR_RETURN (-1);
          }
          src_parts |= FMT_T2P(fmt,rt);
          dst_parts |= FMT_T2P(fmt,ct);
        }

  if (mode==STRICT)
  {
    /* now check whether all parts are covered */
    src_parts &= dst_parts;
    n = BVPD_NPARTS(MG_BVPD(VD_MG(md)));
    for (j=0; j<n; j++)
      if (!(src_parts & (1<<j)))
        REP_ERR_RETURN (-2);
  }
  else if (mode!=NON_STRICT)
    REP_ERR_RETURN (1);

  return (off);
}

/****************************************************************************/
/*D
   MD_nr_nc_mcmpptr_of_ro_co - return comp ptr for row/col object if unique, NULL else

   SYNOPSIS:
   SHORT *MD_nr_nc_mcmpptr_of_ro_co (const MATDATA_DESC *md, INT rowobj, INT colobj, INT *nr, INT *nc, INT mode)

   PARAMETERS:
   .  md - data decsriptor
   .  rowobj - row object type
   .  colobj - col object type
   .  nr - number of rows (may be NULL)
   .  nc - number of cols (may be NULL)
   .  mode - STRICT or NON_STRICT

   DESCRIPTION:
   This function checks whether all components described in 'rowobj'/'colobj'
   are the same for all mtypes using objects of 'rowobj'/'colobj' and returns a component pointer.
   If the components are not unique a NULL is returned.
   If mode is STRICT vectors are required in ALL vobjects of otype.

   CAUTION: it may happen that in parts of the domain vectors in objects of 'otype'
   are not defined at all!

   RETURN VALUE:
   SHORT *
   .n      number of components in objects of 'otype'
   .n      NULL if not unique
   .n      NULL if not the whole domain is covered
   D*/
/****************************************************************************/

SHORT *MD_nr_nc_mcmpptr_of_ro_co_mod (const MATDATA_DESC *md, INT rowobj, INT colobj, INT *nr, INT *nc, INT mode)
{
  FORMAT *fmt;
  SHORT *cptr;
  INT rt,ct,off,rot,cot,nrow,ncol,ncmp,src_parts,dst_parts,i,j,n;

  if (nr!=NULL) *nr = -1;
  if (nc!=NULL) *nc = -1;

  fmt = MGFORMAT(MD_MG(md));
  rot = 1<<rowobj;
  cot = 1<<colobj;
  off = nrow = ncol = src_parts = dst_parts = 0;
  cptr = NULL;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
        if ((rot & FMT_T2O(fmt,rt)) && (cot & FMT_T2O(fmt,ct)))
        {
          if (nrow==0)
          {
            nrow = MD_ROWS_IN_RT_CT(md,rt,ct);
            ncol = MD_COLS_IN_RT_CT(md,rt,ct);
            ncmp = nrow*ncol;
            cptr = MD_MCMPPTR_OF_RT_CT(md,rt,ct);
          }
          else
          {
            if (MD_ROWS_IN_RT_CT(md,rt,ct)!=nrow)
              REP_ERR_RETURN (NULL);
            if (MD_COLS_IN_RT_CT(md,rt,ct)!=ncol)
              REP_ERR_RETURN (NULL);
            for (i=0; i<ncmp; i++)
              if (MD_MCMP_OF_RT_CT(md,rt,ct,i)!=cptr[i])
                REP_ERR_RETURN (NULL);
          }
          src_parts |= FMT_T2P(fmt,rt);
          dst_parts |= FMT_T2P(fmt,ct);
        }

  if (mode==STRICT)
  {
    /* now check whether all parts are covered */
    src_parts &= dst_parts;
    n = BVPD_NPARTS(MG_BVPD(VD_MG(md)));
    for (j=0; j<n; j++)
      if (!(src_parts & (1<<j)))
        REP_ERR_RETURN (NULL);
  }
  else if (mode!=NON_STRICT)
    REP_ERR_RETURN (NULL);

  if (nr!=NULL) *nr = nrow;
  if (nc!=NULL) *nc = ncol;

  return (cptr);
}

/****************************************************************************/
/*D
   MDusesVOTypeOnly - check whether only one vector object type is used

   SYNOPSIS:
   INT MDusesVOTypeOnly (const MATDATA_DESC *md, INT votype)

   PARAMETERS:
   .  md - data decsriptor
   .  votype - vector object type

   DESCRIPTION:
   This function checks whether only one vector object type is used by the
   MATDATA_DESC (root and dest type).

   CAUTION: it may happen that in parts of the domain vectors in objects of 'votype'
   are not defined at all!

   RETURN VALUE:
   SHORT *
   .n      number of components in objects of 'otype'
   .n      NULL if not unique
   D*/
/****************************************************************************/

INT MDusesVOTypeOnly (const MATDATA_DESC *md, INT votype)
{
  FORMAT *fmt;
  INT rt,ct,otp;

  fmt = MGFORMAT(VD_MG(md));
  otp = 1<<votype;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
      if (MD_ISDEF_IN_RT_CT(md,rt,ct))
      {
        if (otp!=FMT_T2O(fmt,rt))
          return(NO);
        if (otp!=FMT_T2O(fmt,ct))
          return(NO);
      }
  return (YES);
}

/****************************************************************************/
/*D
   SwapPartInterfaceData - swap data at domain part interface

   SYNOPSIS:
   INT SwapPartInterfaceData (INT fl, INT tl, SPID_DESC *spid, INT direction)

   PARAMETERS:
   .  fl - from level
   .  tl - to   level
   .  spid - data structure holding the XXXDATA_DESCs to be swapped
   .  direction - forth or back

   DESCRIPTION:
   Suppose the domain is split into two parts coinciding with two subdomains and
   a third part on the interface (inner boundary) between the subdomains.
   Say we want to have only nodal degrees of freedom (dof) all over the domain.
   Say further we want to solve a PDE X with n unknowns in part 0 and a second one (Y)
   with m unknowns in part 1. Then we need n+m dofs in the interface part (number 2)
   because we have to implement boundary equations for both PDEs there.

   Now for actions to be performed on the data of PDE X we may wish to have all
   data with the same offset in node vectors. This is the purpose of
   'SwapPartInterfaceData'.

   The swapped data are described by one or more pairs of VEC/MATDATA_DESCs. The first
   (SPID_VD) of each pair describes the data of PDE X including the interface (and
   thus covers parts 0 and 2 in our example). The second one (SPID_VDI) describes
   the data of PDE X on the interface only (part 2 in our example).
   'SwapPartInterfaceData' with direction SPID_FORTH then swaps the interface data
   (n for a vector in the example) in a way that they occupy the same positions in
   the NODEVEC as the corresponding data do in NODEVECs of part 0.

   After swapping it is possible to perform an action on the data and using the same
   offset in both part 0 and part 2.

   NOTE:
   The SPID_VD descriptors are changed such that they describe the new situation.

   IMPORTANT:
   To restore the original state (this includes reverting the SPID_VD descriptors)
   call 'SwapPartInterfaceData' again (now with SPID_BACK)!

   RETURN VALUE:
   INT
   .n      0: ok
   .n      n: else
   D*/
/****************************************************************************/

INT SwapPartInterfaceData (INT fl, INT tl, SPID_DESC *spid, INT direction)
{
  MULTIGRID *mg;
  VECTOR *vec;
  MATRIX *mat;
  VECDATA_DESC *vd,*vdi;
  MATDATA_DESC *md,*mdi;
  DOUBLE tmp;
  INT tp,tpn,rt,ct,i,j,n,nn,ni,kn,ki,lev,consider_mat;
  INT min,max,mintp;
  INT nv,nm;
  SHORT *vcmp,*mcmp,*pn,*pi;
  SHORT nvn[NVECTYPES],nvi[NVECTYPES];
  SHORT vn[SPID_NVD_MAX*MAX_VEC_COMP],vi[SPID_NVD_MAX*MAX_VEC_COMP];
  SHORT vnoffset[NVECOFFSETS],vioffset[NVECOFFSETS];
  SHORT nmn[NMATTYPES],nmi[NMATTYPES];
  SHORT mn[SPID_NMD_MAX*MAX_MAT_COMP],mi[SPID_NMD_MAX*MAX_MAT_COMP];
  SHORT mnoffset[NMATOFFSETS],mioffset[NMATOFFSETS];
  static INT old_direction=SPID_BACK;

  /* check that direction has changed */
  if (direction==old_direction)
    REP_ERR_RETURN(1);
  old_direction = direction;

  fl = MAX(0,fl);

  if (SPID_NVD(spid)>0)
    mg = VD_MG(SPID_VD(spid,0));
  else if (SPID_NMD(spid)>0)
    mg = MD_MG(SPID_MD(spid,0));
  else
    /* neither vds nor mds defined */
    REP_ERR_RETURN (1);
  ASSERT(fl>=BOTTOMLEVEL(mg));
  ASSERT(tl<=TOPLEVEL(mg));

  /* get interface components to copy from
     and non-interface components to copy to */


  /*****************/
  /* VECDATA_DESCs */
  /*****************/

  tpn = -1;
  nn = ni = 0;
  for (tp=0; tp<NVECTYPES; tp++)
  {
    kn = ki = 0;
    for (i=0; i<SPID_NVD(spid); i++)
      if (VD_ISDEF_IN_TYPE(SPID_VDI(spid,i),tp))
      {
        /* this is the interface */
        vd = SPID_VDI(spid,i);
        n  = VD_NCMPS_IN_TYPE(vd,tp);
        for (j=0; j<n; j++)
        {
          ASSERT(ni<SPID_NVD_MAX*MAX_VEC_COMP);
          vi[ni++] = VD_CMP_OF_TYPE(vd,tp,j);
          ki++;
        }
      }
      else if (VD_ISDEF_IN_TYPE(SPID_VD(spid,i),tp))
      {
        /* this is non-interface */
        vd = SPID_VD(spid,i);
        n  = VD_NCMPS_IN_TYPE(vd,tp);
        for (j=0; j<n; j++)
        {
          ASSERT(nn<SPID_NVD_MAX*MAX_VEC_COMP);
          vn[nn++] = VD_CMP_OF_TYPE(vd,tp,j);
          kn++;
        }
      }
    nvi[tp] = ki;
    nvn[tp] = kn;
  }

  if (nn==0)
  {
    /* there is interface only: take interface type with minimal offset components */
    min = MAX_I;
    for (tp=0; tp<NVECTYPES; tp++)
    {
      max = 0;
      for (i=0; i<SPID_NVD(spid); i++)
        if (VD_ISDEF_IN_TYPE(SPID_VDI(spid,i),tp))
        {
          /* this is the interface */
          vd = SPID_VDI(spid,i);
          n  = VD_NCMPS_IN_TYPE(vd,tp);
          for (j=0; j<n; j++)
            max = MAX(max,VD_CMP_OF_TYPE(vd,tp,j));
        }
      if (max<min)
      {
        mintp = tp;
        min = max;
      }
    }
    tp = mintp;
    kn = 0;
    for (i=0; i<SPID_NVD(spid); i++)
      if (VD_ISDEF_IN_TYPE(SPID_VDI(spid,i),tp))
      {
        /* this is non-interface */
        vd = SPID_VDI(spid,i);
        n  = VD_NCMPS_IN_TYPE(vd,tp);
        for (j=0; j<n; j++)
        {
          ASSERT(nn<SPID_NVD_MAX*MAX_VEC_COMP);
          vn[nn++] = VD_CMP_OF_TYPE(vd,tp,j);
          kn++;
        }
      }
    nvn[tp] = kn;
  }

  ConstructVecOffsets(nvi,vioffset);
  ConstructVecOffsets(nvn,vnoffset);

  /* check that only one non-interface type exists (for uniqueness) */
  vcmp = NULL;
  for (tp=0; tp<NVECTYPES; tp++)
    if (nvn[tp])
    {
      if (vcmp!=NULL)
        REP_ERR_RETURN(1);
      vcmp = vn+vnoffset[tp];
      nv   = nvn[tp];
      tpn  = tp;
    }

  /* check that number of interface and non-interface comps coincide */
  for (tp=0; tp<NVECTYPES; tp++)
    if (nvi[tp]>0)
      if (nvi[tp]!=nv)
        REP_ERR_RETURN(1);

  /* finally change comps of VD on interface to describe data after swapping */
  for (tp=0; tp<NVECTYPES; tp++)
    for (i=0; i<SPID_NVD(spid); i++)
      if (VD_ISDEF_IN_TYPE(SPID_VDI(spid,i),tp))
      {
        /* this is the interface */
        vdi = SPID_VDI(spid,i);
        vd  = SPID_VD(spid,i);
        n   = VD_NCMPS_IN_TYPE(vd,tp);
        ASSERT(n==VD_NCMPS_IN_TYPE(vd,tpn));
        if (direction==SPID_FORTH)
          /* get comps from non-interface */
          for (j=0; j<n; j++)
            VD_CMP_OF_TYPE(vd,tp,j) = VD_CMP_OF_TYPE(vd,tpn,j);
        else if (direction==SPID_BACK)
          /* get comps from interface */
          for (j=0; j<n; j++)
            VD_CMP_OF_TYPE(vd,tp,j) = VD_CMP_OF_TYPE(vdi,tp,j);
        else
          REP_ERR_RETURN(1);
      }

  /*****************/
  /* MATDATA_DESCs */
  /*****************/

  nn = ni = 0;
  for (tp=0; tp<NMATTYPES; tp++)
  {
    kn = ki = 0;
    for (i=0; i<SPID_NMD(spid); i++)
      if (MD_ISDEF_IN_MTYPE(SPID_MDI(spid,i),tp))
      {
        /* this is the interface */
        md = SPID_MDI(spid,i);
        n  = MD_NCMPS_IN_MTYPE(md,tp);
        for (j=0; j<n; j++)
        {
          ASSERT(ni<SPID_NMD_MAX*MAX_MAT_COMP);
          mi[ni++] = MD_MCMP_OF_MTYPE(md,tp,j);
          ki++;
        }
      }
      else if (MD_ISDEF_IN_MTYPE(SPID_MD(spid,i),tp))
      {
        /* this is non-interface */
        md = SPID_MD(spid,i);
        n  = MD_NCMPS_IN_MTYPE(md,tp);
        for (j=0; j<n; j++)
        {
          ASSERT(nn<SPID_NMD_MAX*MAX_MAT_COMP);
          mn[nn++] = MD_MCMP_OF_MTYPE(md,tp,j);
          kn++;
        }
      }
    nmi[tp] = ki;
    nmn[tp] = kn;
  }

  if (nn==0)
  {
    /* there is interface only: take interface type with minimal offset components */
    min = MAX_I;
    for (tp=0; tp<NMATTYPES; tp++)
    {
      max = 0;
      for (i=0; i<SPID_NMD(spid); i++)
        if (MD_ISDEF_IN_MTYPE(SPID_MDI(spid,i),tp))
        {
          /* this is the interface */
          md = SPID_MDI(spid,i);
          n  = MD_NCMPS_IN_MTYPE(md,tp);
          for (j=0; j<n; j++)
            max = MAX(max,MD_MCMP_OF_MTYPE(md,tp,j));
        }
      if (max<min)
      {
        mintp = tp;
        min = max;
      }
    }
    tp = mintp;
    kn = 0;
    for (i=0; i<SPID_NMD(spid); i++)
      if (MD_ISDEF_IN_MTYPE(SPID_MDI(spid,i),tp))
      {
        /* this is non-interface */
        md = SPID_MDI(spid,i);
        n  = MD_NCMPS_IN_MTYPE(md,tp);
        for (j=0; j<n; j++)
        {
          ASSERT(nn<SPID_NMD_MAX*MAX_MAT_COMP);
          mn[nn++] = MD_MCMP_OF_MTYPE(md,tp,j);
          kn++;
        }
      }
    nmn[tp] = kn;
  }

  ConstructMatOffsetsAlt(nmi,mioffset);
  ConstructMatOffsetsAlt(nmn,mnoffset);

  /* check that only one non-interface type exists (for uniqueness) */
  mcmp = NULL;
  for (tp=0; tp<NMATTYPES; tp++)
    if (nmn[tp])
    {
      ASSERT(mcmp==NULL);
      mcmp = mn+mnoffset[tp];
      nm   = nmn[tp];
      if (tpn!=-1)
        ASSERT(tp==tpn*NVECTYPES+tpn);
      tpn  = tp;
    }

  /* check that number of interface and non-interface comps coincide */
  for (tp=0; tp<NMATTYPES; tp++)
    if (nmi[tp]>0)
      if (nmi[tp]!=nm)
        REP_ERR_RETURN(1);

  /* finally change comps of MD on interface to describe data after swapping */
  for (tp=0; tp<NMATTYPES; tp++)
    for (i=0; i<SPID_NMD(spid); i++)
      if (MD_ISDEF_IN_MTYPE(SPID_MDI(spid,i),tp))
      {
        /* this is the interface */
        mdi = SPID_MDI(spid,i);
        md  = SPID_MD(spid,i);
        n   = MD_NCMPS_IN_MTYPE(md,tp);
        ASSERT(n==MD_NCMPS_IN_MTYPE(md,tpn));
        if (direction==SPID_FORTH)
          /* get comps from non-interface */
          for (j=0; j<n; j++)
            MD_MCMP_OF_MTYPE(md,tp,j) = MD_MCMP_OF_MTYPE(md,tpn,j);
        else if (direction==SPID_BACK)
          /* get comps from interface */
          for (j=0; j<n; j++)
            MD_MCMP_OF_MTYPE(md,tp,j) = MD_MCMP_OF_MTYPE(mdi,tp,j);
        else
          REP_ERR_RETURN(1);
      }


  /*********************************/
  /* now loop vectors and matrices */
  /*              and swap data			 */
  /*********************************/

  consider_mat = (SPID_NMD(spid)>0);
  if (direction==SPID_FORTH)
    for (lev=fl; lev<=tl; lev++)
      for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!=NULL; vec=SUCCVC(vec))
      {
        rt = VTYPE(vec);
        n  = nvi[rt];
        if (n>0)
        {
          pn = vcmp;
          pi = vi+vioffset[rt];
          for (i=0; i<nv; i++,pi++,pn++)
            SWAP_VEC_DATA(vec,pi,pn);
        }
        if (consider_mat)
          for (mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat))
          {
            ct = VTYPE(MDEST(mat));
            tp = MTP(rt,ct);
            n  = nmi[tp];
            if (n>0)
            {
              pn = mcmp;
              pi = mi+mioffset[tp];
              for (i=0; i<nm; i++,pi++,pn++)
                SWAP_MAT_DATA(mat,pi,pn);
            }
          }
      }
  else if (direction==SPID_BACK)
    /* run swapping loops backwards */
    for (lev=fl; lev<=tl; lev++)
      for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!=NULL; vec=SUCCVC(vec))
      {
        rt = VTYPE(vec);
        n  = nvi[rt];
        if (n>0)
        {
          pn = vcmp+nv-1;
          pi = vi+vioffset[rt]+nv-1;
          for (i=nv-1; i>=0; i--,pi--,pn--)
            SWAP_VEC_DATA(vec,pi,pn);
        }
        if (consider_mat)
          for (mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat))
          {
            ct = VTYPE(MDEST(mat));
            tp = MTP(rt,ct);
            n  = nmi[tp];
            if (n>0)
            {
              pn = mcmp+nm-1;
              pi = mi+mioffset[tp]+nm-1;
              for (i=nm-1; i>=0; i--,pi--,pn--)
                SWAP_MAT_DATA(mat,pi,pn);
            }
          }
      }
  else
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/*D
   SwapPartSkipflags - swap the skip bits of vectors

   SYNOPSIS:
   INT SwapPartSkipflags (INT fl, INT tl, const VECDATA_DESC *vdg, const VECDATA_DESC *vdi, INT direction)

   PARAMETERS:
   .  fl - from level
   .  tl - to   level
   .  vdg - descriptor for global data
   .  vdi - descriptor for part data restricted to interface
   .  direction - forth or back

   DESCRIPTION:
   Similar to 'SwapPartInterfaceData' this function swaps the skip bits of vectors
   such that the bits of the part descriptor ly at the lowest bits of the vectors
   skip component.

   NOTE:
   Even if several data descriptors are swapped by 'SwapPartInterfaceData' swap skip bits
   only once!

   IMPORTANT:
   To restore the original state call 'SwapPartSkipflags' again (now with SPID_BACK)!

   RETURN VALUE:
   INT
   .n      0: ok
   .n      n: else
   D*/
/****************************************************************************/

INT SwapPartSkipflags (INT fl, INT tl, const VECDATA_DESC *vdg, const VECDATA_DESC *vdi, INT direction)
{
  MULTIGRID *mg;
  VECTOR *vec;
  struct {
    INT len;                                            /* number of bits		*/
    INT shift;                                          /* shift in global desc	*/
    unsigned INT bits;                          /* mask for bits		*/
    unsigned INT rest;                          /* mask for rest		*/
  } sps[NVECTYPES],*p;
  INT i,j,n,tp,shift,cmpi,lev,skip,bits,rest;

  mg = VD_MG(vdg);

  fl = MAX(0,fl);

  /* first compute mask, len and offset of vdi skipbits on interface */
  for (tp=0; tp<NVECTYPES; tp++)
  {
    sps[tp].len = 0;
    if (VD_ISDEF_IN_TYPE(vdi,tp))
    {
      ASSERT(VD_ISDEF_IN_TYPE(vdg,tp));

      n = VD_NCMPS_IN_TYPE(vdi,tp);

      /* compute position of first vdi comp rel to vdg */
      cmpi = VD_CMP_OF_TYPE(vdi,tp,0);
      for (j=0; j<VD_NCMPS_IN_TYPE(vdg,tp); j++)
        if (VD_CMP_OF_TYPE(vdg,tp,j)==cmpi)
          break;
      ASSERT(j<VD_NCMPS_IN_TYPE(vdg,tp));
      shift = j;

      /* are the components of vdi subsequent in vdg? */
      if (shift+n>VD_NCMPS_IN_TYPE(vdg,tp))
        REP_ERR_RETURN(1);
      for (i=1; i<n; i++)
        if (VD_CMP_OF_TYPE(vdi,tp,i)!=VD_CMP_OF_TYPE(vdg,tp,shift+i))
          REP_ERR_RETURN(1);

      if (shift==0)
        /* there will be nothing to swap */
        continue;

      /* fill sps */
      sps[tp].len             = n;
      sps[tp].shift   = shift;

      sps[tp].bits    = POW2(n)-1;
      if (direction==SPID_FORTH)
        sps[tp].bits    = sps[tp].bits << shift;

      sps[tp].rest    = ~(sps[tp].bits);
    }
  }

  /* now loop vectors */
  if (direction==SPID_FORTH)
  {
    for (lev=fl; lev<=tl; lev++)
      for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!=NULL; vec=SUCCVC(vec))
        if (sps[VTYPE(vec)].len)
        {
          p = sps+VTYPE(vec);

          skip = VECSKIP(vec);
          if (skip==0)
            continue;
          bits = skip & p->bits;
          rest = skip & p->rest;
          VECSKIP(vec) = (bits >> p->shift) | (rest << p->shift);
        }
  }
  else if (direction==SPID_BACK)
  {
    for (lev=fl; lev<=tl; lev++)
      for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!=NULL; vec=SUCCVC(vec))
        if (sps[VTYPE(vec)].len)
        {
          p = sps+VTYPE(vec);

          skip = VECSKIP(vec);
          if (skip==0)
            continue;
          bits = skip & p->bits;
          rest = skip & p->rest;
          VECSKIP(vec) = (bits << p->shift) | (rest >> p->shift);
        }
  }
  else
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/*
   InitUserDataManager - Init this file

   SYNOPSIS:
   INT InitUserDataManager (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    __LINE__ if error occured.
 */
/****************************************************************************/

INT InitUserDataManager ()
{
  char *names;
  INT i;

  MatrixDirID = GetNewEnvDirID();
  VectorDirID = GetNewEnvDirID();
  MatrixVarID = GetNewEnvVarID();
  VectorVarID = GetNewEnvVarID();

  names = DEFAULT_NAMES;

  for (i=0; i<MAX(MAX_VEC_COMP,strlen(DEFAULT_NAMES)); i++)
    NoVecNames[i] = names[i];
  for (i=0; i<2*MAX_MAT_COMP; i++)
    NoMatNames[i] = ' ';

  return (NUM_OK);
}
