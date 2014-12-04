// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      transgrid.c                                                   */
/*                                                                          */
/* Purpose:   standard grid transfer functions (restriction/interpolation)  */
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   20.04.95 begin, ug version 3.0								*/
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

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "ugtypes.h"
#include "ugdevices.h"
#include "misc.h"
#include "gm.h"
#include "elements.h"
#include "algebra.h"
#include "ugdevices.h"
#include "evm.h"
#include "shapes.h"
#include "debug.h"
#include "general.h"
#include "block.h"
#include "cw.h"

#include "np.h"
#include "disctools.h"
#include "transgrid.h"

USING_UG_NAMESPACES

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

enum {ZERO_BLOCK, DIAG_BLOCK, FULL_BLOCK};

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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/** \brief Restrict defect of fine node vectors

   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor
   \param damp - damping factor for every component

   This function restricts defect of fine node vectors with NEWDEFECT_CLASS
   to the next coarser grid. It is the transposent operation to
   'StandardIntCorNodeVector'.
   First, it resets all components to zero. It considers the VECSKIP-flags.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

static INT StandardRestrictNodeVector (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp)
{
  GRID *CoarseGrid;
  ELEMENT *theElement;
  VERTEX *theVertex;
  NODE *theNode;
  EDGE *theEdge;
  VECTOR *v,*vc;
  DOUBLE c[MAX_CORNERS_OF_ELEM],s[MAX_SINGLE_VEC_COMP];
  const SHORT *toComp,*fromComp,*edComp;
  INT i,j,n,ncomp,edcomp,vecskip,dt;

  CoarseGrid = DOWNGRID(FineGrid);

  toComp    = VD_ncmp_cmpptr_of_otype_mod(to,NODEVEC,&ncomp,NON_STRICT);
  fromComp  = VD_cmpptr_of_otype_mod(from,NODEVEC,NON_STRICT);
  edComp   = VD_ncmp_cmpptr_of_otype_mod(to,EDGEVEC,&edcomp,NON_STRICT);
  if (ncomp <= 0)
    return(NUM_ERROR);
  if (ncomp>MAX_SINGLE_VEC_COMP)
    return (NUM_BLOCK_TOO_LARGE);
  if (ncomp < edcomp)
    return(NUM_ERROR);

  /* reset coarser defect at positions where a new defect is restricted */
  dt = VD_DATA_TYPES(to);
  for (v=PFIRSTVECTOR(CoarseGrid); v!= NULL; v=SUCCVC(v))
    if (V_IN_DATATYPE(v,dt))
      if (VNCLASS(v)>=NEWDEF_CLASS)
        for (i=0; i<ncomp; i++)
          VVALUE(v,toComp[i]) = 0.0;

  /* compute contributions to all coarse node vectors */
  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    v = NVECTOR(theNode);
    if (VCLASS(v)<NEWDEF_CLASS) continue;
    if (!V_IN_DATATYPE(v,dt)) continue;
#ifdef __PERIODIC_BOUNDARY__
    if ((NODE *)VOBJECT(v)!=theNode) continue;
#endif
    if (CORNERTYPE(theNode))
    {
      vc = NVECTOR((NODE *)NFATHER(theNode));
      vecskip  = VECSKIP(vc);
      for (i=0; i<ncomp; i++)
        if (!(vecskip & (1<<i)))
          VVALUE(vc,toComp[i]) += damp[i] * VVALUE(v,fromComp[i]);
    }
    else
    {
      theVertex  = MYVERTEX(theNode);
      theElement = VFATHER(theVertex);
      n = CORNERS_OF_ELEM(theElement);
      GNs(n,LCVECT(theVertex),c);
      for (i=0; i<ncomp; i++)
        s[i] = damp[i] * VVALUE(v,fromComp[i]);
      for (i=0; i<n; i++)
      {
        vc = NVECTOR(CORNER(theElement,i));

        /* HRR_TODO: at coupling boundaries we want to avoid accessing vectors of
           other parts. This may introduce a small error at curved inner coupling
           boundaries where the vertex does not ly on the element edge */
        if (!V_IN_DATATYPE(vc,dt)) continue;

        vecskip = VECSKIP(vc);
        for (j=0; j<ncomp; j++)
          if (!(vecskip & (1<<j)))
            VVALUE(vc,toComp[j]) += c[i] * s[j];
      }
    }
  }

  for (v=PFIRSTVECTOR(FineGrid); v!= NULL; v=SUCCVC(v)) {
    if (VOTYPE(v) != EDGEVEC) continue;
    theEdge = (EDGE *)VOBJECT(v);
    theNode = NBNODE(LINK0(theEdge));
    if (NTYPE(theNode) == MID_NODE)
      vc = EDVECTOR(NFATHEREDGE(theNode));
    else if (NTYPE(theNode) == CORNER_NODE)
      vc = NVECTOR(NFATHER(theNode));
    vecskip = VECSKIP(vc);
    for (j=0; j<edcomp; j++)
      if (!(vecskip & (1<<j)))
        VVALUE(vc,toComp[j]) += 0.5 * damp[j] * VVALUE(v,edComp[j]);
    theNode = NBNODE(LINK1(theEdge));
    if (NTYPE(theNode) == MID_NODE)
      vc = EDVECTOR(NFATHEREDGE(theNode));
    else if (NTYPE(theNode) == CORNER_NODE)
      vc = NVECTOR(NFATHER(theNode));
    vecskip = VECSKIP(vc);
    for (j=0; j<edcomp; j++)
      if (!(vecskip & (1<<j)))
        VVALUE(vc,toComp[j]) += 0.5 * damp[j] * VVALUE(v,edComp[j]);
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Interpolate correction from coarse node vectors

   \param FineGrid - pointer to grid
   \param to - vector descriptor
   \param from  - vector descriptor
   \param damp - damping factor for every component

   This function interpolates correction from coarse node vectors,
   using linear resp. bilinear interpolation.
   First, it resets all components to zero. It considers the VECSKIP-flags.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

static INT StandardIntCorNodeVector (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp)
{
  GRID *CoarseGrid;
  ELEMENT *theElement;
  VERTEX *theVertex;
  NODE *theNode;
  EDGE *theEdge;
  VECTOR *v,*vc,*cvec[MAX_CORNERS_OF_ELEM];
  DOUBLE c[MAX_CORNERS_OF_ELEM];
  const SHORT *toComp,*fromComp,*edComp;
  INT i,j,n,ncomp,edcomp,vecskip,skip,dt;

  CoarseGrid = DOWNGRID(FineGrid);

  toComp   = VD_ncmp_cmpptr_of_otype_mod(to,NODEVEC,&ncomp,NON_STRICT);
  edComp   = VD_ncmp_cmpptr_of_otype_mod(to,EDGEVEC,&edcomp,NON_STRICT);
  fromComp = VD_cmpptr_of_otype_mod(from,NODEVEC,NON_STRICT);
  if (ncomp <= 0)
    return(NUM_ERROR);
  if (ncomp < edcomp)
    return(NUM_ERROR);

  /* reset fine to field */
  dt = VD_DATA_TYPES(to);
  for (v=FIRSTVECTOR(FineGrid); v!= NULL; v=SUCCVC(v))
    if (V_IN_DATATYPE(v,dt))
      for (i=0; i<ncomp; i++)
        VVALUE(v,toComp[i]) = 0.0;

  /* compute contributions from all coarse node vectors */
  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    v = NVECTOR(theNode);
    if (!V_IN_DATATYPE(v,dt)) continue;
#ifdef __PERIODIC_BOUNDARY__
    if ((NODE *)VOBJECT(v)!=theNode) continue;
#endif
    vecskip = VECSKIP(v);
    skip = TRUE;
    for (i=0; i<ncomp; i++)
      if (!(vecskip & (1<<i)))
        skip = FALSE;
    if (skip) continue;                         /* skip only if all flags are TRUE */

    if (CORNERTYPE(theNode))
    {
      vc = NVECTOR((NODE *)NFATHER(theNode));
      for (i=0; i<ncomp; i++)
        if (!(vecskip & (1<<i)))
          VVALUE(v,toComp[i]) = damp[i] * VVALUE(vc,fromComp[i]);
    }
    else
    {
      theVertex  = MYVERTEX(theNode);
      theElement = VFATHER(theVertex);
      n = CORNERS_OF_ELEM(theElement);
      GNs(n,LCVECT(theVertex),c);
      for (i=0; i<n; i++)
      {
        cvec[i] = NVECTOR(CORNER(theElement,i));
        /* HRR_TODO: at coupling boundaries we want to avoid accessing vectors of
           other parts. This may introduce a small error at curved inner coupling
           boundaries where the vertex does not ly on the element edge */
        if (!V_IN_DATATYPE(cvec[i],dt))
          c[i] = 0.0;
      }
      for (j=0; j<ncomp; j++)
        if (!(vecskip & (1<<j)))
          for (i=0; i<n; i++)
            if (c[i]!=0.0)                                      /* since VVALUE(cvec[i],fromComp[j]) maybe not valid */
              VVALUE(v,toComp[j]) += c[i]*damp[j] * VVALUE(cvec[i],fromComp[j]);
    }
  }

  if (edcomp <= 0)
    return (NUM_OK);

  for (v=PFIRSTVECTOR(FineGrid); v!= NULL; v=SUCCVC(v)) {
    if (VOTYPE(v) != EDGEVEC) continue;
    theEdge = (EDGE *)VOBJECT(v);
    theNode = NBNODE(LINK0(theEdge));
    if (NTYPE(theNode) == MID_NODE)
      vc = EDVECTOR(NFATHEREDGE(theNode));
    else if (NTYPE(theNode) == CORNER_NODE)
      vc = NVECTOR(NFATHER(theNode));
    vecskip = VECSKIP(vc);
    for (j=0; j<edcomp; j++)
      if (!(vecskip & (1<<j)))
        VVALUE(v,edComp[j]) += 0.5 * damp[j] * VVALUE(vc,fromComp[j]);
    theNode = NBNODE(LINK1(theEdge));
    if (NTYPE(theNode) == MID_NODE)
      vc = EDVECTOR(NFATHEREDGE(theNode));
    else if (NTYPE(theNode) == CORNER_NODE)
      vc = NVECTOR(NFATHER(theNode));
    vecskip = VECSKIP(vc);
    for (j=0; j<edcomp; j++)
      if (!(vecskip & (1<<j)))
        VVALUE(v,edComp[j]) += 0.5 * damp[j] * VVALUE(vc,fromComp[j]);
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Interpolate the solution to the new vectors

   \param FineGrid - pointer to grid
   \param Cor - vector descriptor

   This function interpolates the solution from coarse node vectors
   to new vectors, using linear resp. bilinear interpolation and
   considering the VECSKIP-flags.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

static INT StandardIntNewNodeVector (GRID *FineGrid, const VECDATA_DESC *Cor)
{
  GRID *CoarseGrid;
  ELEMENT *theElement;
  VERTEX *theVertex;
  NODE *theNode;
  EDGE *theEdge;
  VECTOR *v,*vc,*cvec[MAX_CORNERS_OF_ELEM];
  DOUBLE c[MAX_CORNERS_OF_ELEM];
  const SHORT *Comp,*edComp;
  INT i,j,n,ncomp,edcomp,dt;

  CoarseGrid = DOWNGRID(FineGrid);

  Comp   = VD_ncmp_cmpptr_of_otype_mod(Cor,NODEVEC,&ncomp,NON_STRICT);
  edComp   = VD_ncmp_cmpptr_of_otype_mod(Cor,EDGEVEC,&edcomp,NON_STRICT);
  if (ncomp <= 0)
    return(NUM_ERROR);
  if (ncomp < edcomp)
    return(NUM_ERROR);

  /* interpolate values to all fine node vectors */
  dt = VD_DATA_TYPES(Cor);
  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    v = NVECTOR(theNode);
#ifdef __PERIODIC_BOUNDARY__
    if ((NODE *)VOBJECT(v)!=theNode) continue;
#endif
    if (!VNEW(v)) continue;
    if (!V_IN_DATATYPE(v,dt)) continue;

    if (CORNERTYPE(theNode))
    {
      vc = NVECTOR((NODE *)NFATHER(theNode));
      for (i=0; i<ncomp; i++)
        VVALUE(v,Comp[i]) = VVALUE(vc,Comp[i]);
    }
    else
    {
      theVertex  = MYVERTEX(theNode);
      theElement = VFATHER(theVertex);
      n   = CORNERS_OF_ELEM(theElement);
      GNs(n,LCVECT(theVertex),c);
      for (i=0; i<n; i++)
      {
        cvec[i] = NVECTOR(CORNER(theElement,i));
        /* HRR_TODO: at coupling boundaries we want to avoid accessing vectors of
           other parts. This may introduce a small error at curved inner coupling
           boundaries where the vertex does not lie on the element edge */
        if (!V_IN_DATATYPE(cvec[i],dt))
          c[i] = 0.0;
      }
      for (j=0; j<ncomp; j++)
      {
        VVALUE(v,Comp[j]) = 0.0;

        for (i=0; i<n; i++)
          if (c[i]!=0.0)                                /* since VVALUE(cvec[i],fromComp[j]) maybe not valid */
            VVALUE(v,Comp[j]) += c[i] * VVALUE(cvec[i],Comp[j]);
      }
    }
  }

  if (edcomp <= 0)
    return (NUM_OK);

  for (v=PFIRSTVECTOR(FineGrid); v!= NULL; v=SUCCVC(v)) {
    if (!VNEW(v)) continue;
    if (VOTYPE(v) != EDGEVEC) continue;
    theEdge = (EDGE *)VOBJECT(v);
    theNode = NBNODE(LINK0(theEdge));
    if (NTYPE(theNode) == MID_NODE)
      vc = EDVECTOR(NFATHEREDGE(theNode));
    else if (NTYPE(theNode) == CORNER_NODE)
      vc = NVECTOR(NFATHER(theNode));
    for (j=0; j<edcomp; j++)
      VVALUE(v,edComp[j]) += 0.5 * VVALUE(vc,Comp[j]);
    theNode = NBNODE(LINK1(theEdge));
    if (NTYPE(theNode) == MID_NODE)
      vc = EDVECTOR(NFATHEREDGE(theNode));
    else if (NTYPE(theNode) == CORNER_NODE)
      vc = NVECTOR(NFATHER(theNode));
    for (j=0; j<edcomp; j++)
      VVALUE(v,edComp[j]) += 0.5 * VVALUE(vc,Comp[j]);
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Restrict defect of fine vectors with NEWDEFECT_CLASS

   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor
   \param damp - damping factor for every component

   This function restricts defect of fine vectors with NEWDEFECT_CLASS,
   considers the VECSKIP-flags.
   It calls 'StandardRestrictNodeVector'.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX StandardRestrict (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp)
{
  FORMAT *fmt;
  INT vtype,rv,otype;
  const SHORT *offset;

  if (DOWNGRID(FineGrid)==NULL)
    return (NUM_NO_COARSER_GRID);

  if (0 < VD_NCMPS_IN_TYPE(to,EDGEVEC))
    if (VD_NCMPS_IN_TYPE(to,EDGEVEC) <= VD_NCMPS_IN_TYPE(to,NODEVEC)) {
      rv = StandardRestrictNodeVector(FineGrid,to,from,
                                      (const DOUBLE *)damp);
      return (rv);
    }

  offset = VD_OFFSETPTR(to);
  fmt = MGFORMAT(MYMG(FineGrid));

  for (otype=0; otype<MAXVOBJECTS; otype++)
    if (VD_OBJ_USED(to) & BITWISE_TYPE(otype))
      switch (otype)
      {
      case ELEMVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      case NODEVEC :
        for (vtype=0; vtype<NVECTYPES; vtype++)
          if (VD_ISDEF_IN_TYPE(to,vtype))
            if (GetUniqueOTypeOfVType(fmt,vtype)<0)
              REP_ERR_RETURN(1)
              if ((rv=StandardRestrictNodeVector(FineGrid,to,from,damp+offset[otype]))!=NUM_OK)
                return (rv);
        break;
      case EDGEVEC :
        if ((rv=StandardRestrictNodeVector(FineGrid,to,from,damp+offset[otype]))!=NUM_OK)
          return (rv);
        break;
      case SIDEVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      default :
        UserWrite("not unique");
        return (NUM_ERROR);
      }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Interpolates the correction of the coarse grids

   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor
   \param damp - damping factor for every component

   This function interpolates correction from coarse side vectors,
   considers the VECSKIP-flags.
   It calls 'StandardIntCorNodeVector'.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX StandardInterpolateCorrection (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp)
{
  FORMAT *fmt;
  INT vtype,rv,otype;
  const SHORT *offset;

  if (DOWNGRID(FineGrid)==NULL)
    return (NUM_NO_COARSER_GRID);

  if (0 < VD_NCMPS_IN_TYPE(to,EDGEVEC))
    if (VD_NCMPS_IN_TYPE(to,EDGEVEC) < VD_NCMPS_IN_TYPE(to,NODEVEC)) {
      rv=StandardIntCorNodeVector(FineGrid,to,from,(const DOUBLE *)damp);
      return (rv);
    }

  offset = VD_OFFSETPTR(to);
  fmt = MGFORMAT(MYMG(FineGrid));

  for (otype=0; otype<MAXVOBJECTS; otype++)
    if (VD_OBJ_USED(to) & BITWISE_TYPE(otype))
      switch (otype)
      {
      case ELEMVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      case NODEVEC :
        for (vtype=0; vtype<NVECTYPES; vtype++)
          if (VD_ISDEF_IN_TYPE(to,vtype))
            if (GetUniqueOTypeOfVType(fmt,vtype)<0)
              REP_ERR_RETURN(1)
              if ((rv=StandardIntCorNodeVector(FineGrid,to,from,damp+offset[otype]))!=NUM_OK)
                return (rv);
        break;
      case EDGEVEC :
        if ((rv=StandardIntCorNodeVector(FineGrid,to,from,damp+offset[otype]))!=NUM_OK)
          return (rv);
        break;
      case SIDEVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      default :
        UserWrite("not unique");
        return (NUM_ERROR);
      }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Interpolates the solution on the new vectors

   \param FineGrid - pointer to grid
   \param Sol - type vector descriptor

   This function interpolates the solution from coarse vectors
   to new vectors, considering the VECSKIP-flags.
   It calls 'StandardIntNewNodeVector'.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX StandardInterpolateNewVectors (GRID *FineGrid, const VECDATA_DESC *Sol)
{
  FORMAT *fmt;
  INT vtype,rv,otype;

  if (DOWNGRID(FineGrid)==NULL)
    return (NUM_NO_COARSER_GRID);

  if (0 < VD_NCMPS_IN_TYPE(Sol,EDGEVEC))
    if (VD_NCMPS_IN_TYPE(Sol,EDGEVEC) < VD_NCMPS_IN_TYPE(Sol,NODEVEC)) {
      rv = StandardIntNewNodeVector(FineGrid,Sol);
      return (rv);
    }

  fmt = MGFORMAT(MYMG(FineGrid));

  for (otype=0; otype<MAXVOBJECTS; otype++)
    if (VD_OBJ_USED(Sol) & BITWISE_TYPE(otype))
      switch (otype)
      {
      case ELEMVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      case NODEVEC :
        for (vtype=0; vtype<NVECTYPES; vtype++)
          if (VD_ISDEF_IN_TYPE(Sol,vtype))
            if (GetUniqueOTypeOfVType(fmt,vtype)<0)
              REP_ERR_RETURN(1)
              if ((rv=StandardIntNewNodeVector(FineGrid,Sol))!=NUM_OK)
                return (rv);
        break;
      case EDGEVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      case SIDEVEC :
        UserWrite("not implemented");
        return (NUM_ERROR);
      default :
        UserWrite("not unique");
        return (NUM_ERROR);
      }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Project node values on lower levels

   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor

   This function projects node values to the father node resp. to
   edge values (depending on the format).

   \return
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX StandardProject (GRID *CoarseGrid, const VECDATA_DESC *to,
                                   const VECDATA_DESC *from)
{
  VECTOR *v;
  NODE *theNode;
  DOUBLE *val;
  const SHORT *toComp,*fromComp,*edComp;
  INT i,ncomp,edcomp,nfrom,dt;
    #ifndef ModelP
  ELEMENT *t;
  VECTOR *v0[MAX_NODAL_VECTORS],*v1[MAX_NODAL_VECTORS];
  INT j,m;
        #endif

  toComp   = VD_ncmp_cmpptr_of_otype_mod(to,NODEVEC,&ncomp,NON_STRICT);
  edComp   = VD_ncmp_cmpptr_of_otype_mod(to,EDGEVEC,&edcomp,NON_STRICT);
  fromComp = VD_ncmp_cmpptr_of_otype_mod(from,NODEVEC,&nfrom,NON_STRICT);

  if (ncomp <= 0)
    return (NUM_OK);
  if (ncomp < edcomp)
    return(NUM_ERROR);
  if (nfrom < ncomp)
    return(NUM_ERROR);
  if (ncomp>MAX_VEC_COMP)
    return (NUM_BLOCK_TOO_LARGE);

  dt = VD_DATA_TYPES(to);
  for (v=PFIRSTVECTOR(CoarseGrid); v!= NULL; v=SUCCVC(v))
    if (V_IN_DATATYPE(v,dt))
    {
      if (VOTYPE(v) == NODEVEC)
      {
        theNode = SONNODE((NODE *)VOBJECT(v));
        if (theNode ==  NULL)
          continue;
        val = VVALUEPTR(NVECTOR(theNode),0);
        for (i=0; i<ncomp; i++)
          VVALUE(v,toComp[i]) = val[fromComp[i]];
      }
      else if (VOTYPE(v) == EDGEVEC)
      {
        theNode = MIDNODE((EDGE *)VOBJECT(v));
        if (theNode ==  NULL)
          continue;
        val = VVALUEPTR(NVECTOR(theNode),0);
        for (i=0; i<edcomp; i++)
          VVALUE(v,edComp[i]) = val[fromComp[i]];
      }
    }

  if (edcomp <= 0)
    return (NUM_OK);

    #ifndef ModelP
  fromComp = VD_cmpptr_of_otype(from,EDGEVEC);
  for (t=FIRSTELEMENT(CoarseGrid); t!=NULL; t=SUCCE(t)) {
    if (NSONS(t) != 1)
      continue;
    GetVectorsOfEdges ((const ELEMENT *)t       ,&m,v0);
    DataTypeFilterVList(dt,v0,&m);
    GetVectorsOfEdges ((const ELEMENT *)SON(t,0),&m,v1);
    DataTypeFilterVList(dt,v1,&m);
    for (j=0; j<m; j++)
      for (i=0; i<edcomp; i++)
        VVALUE(v0[j],edComp[i]) = VVALUE(v1[j],fromComp[i]);
  }
    #endif

  return (NUM_OK);
}

#ifdef __INTERPOLATION_MATRIX__

/****************************************************************************/
/** \brief Set all interpolation matrix entries to 0

   \param g - pointer to a grid
   \param theVD - vector descriptor

   This function sets all interpolation matrix entries to 0 and
   sets VINDEX to 0 for all vectors.

   \return
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR 1 if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX ClearIMatrix (GRID *g, VECDATA_DESC *theVD)
{
  VECTOR *v;
  MATRIX *m;

  if (VD_IS_SCALAR(theVD)) {
    register SHORT mask = VD_SCALTYPEMASK(theVD);

    for (v = PFIRSTVECTOR(g); v != NULL; v = SUCCVC(v)) {
      VINDEX(v) = 0;
      if (VDATATYPE(v)&mask) {
        for (m = VISTART(v); m != NULL; m = NEXT(m))
          MVALUE(m,0) = 0.0;
      }
    }
    return (NUM_OK);
  }
  for (v = PFIRSTVECTOR(g); v != NULL; v = SUCCVC(v)) {
    INT j,rcomp;
    DOUBLE *mptr;

    VINDEX(v) = 0;
    rcomp = VD_NCMPS_IN_TYPE(theVD,VTYPE(v));
    for (m = VISTART(v); m != NULL; m = NEXT(m)) {
      mptr = MVALUEPTR(m,0);
      for (j=0; j<rcomp*VD_NCMPS_IN_TYPE(theVD,MDESTTYPE(m)); j++)
        mptr[j] = 0.0;
    }
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Scale the interpolation matrix

   \param g - pointer to a grid
   \param theVD - vector descriptor

   This function scales all interpolation matrix entries by 1 / VINDEX
   resets the VINDEX for all vectors.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR 1 if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX ScaleIMatrix (GRID *g, VECDATA_DESC *theVD)
{
  VECTOR *v;
  MATRIX *m;
  INT i,j,rcomp;
  DOUBLE scale,*mptr;

  i = 0;
  if (VD_IS_SCALAR(theVD))
  {
    for (v = PFIRSTVECTOR(g); v != NULL; v = SUCCVC(v))
    {
      if (VINDEX(v) > 1)
      {
        scale = 1.0 / VINDEX(v);
        for (m = VISTART(v); m != NULL; m = NEXT(m))
          MVALUE(m,0) *= scale;
      }
      VINDEX(v) = i++;
    }
    return (NUM_OK);
  }

  for (v = PFIRSTVECTOR(g); v != NULL; v = SUCCVC(v))
  {
    if (VINDEX(v) > 1)
    {
      scale = 1.0 / VINDEX(v);
      rcomp = VD_NCMPS_IN_TYPE(theVD,VTYPE(v));
      for (m = VISTART(v); m != NULL; m = NEXT(m))
      {
        mptr = MVALUEPTR(m,0);
        for (j=0; j<rcomp*VD_NCMPS_IN_TYPE(theVD,MDESTTYPE(m)); j++)
          mptr[j] *= scale;
      }
    }
    VINDEX(v) = i++;
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Reset vindex

   \param g - pointer to a grid

   This function sets VINDEX to 0 for all vectors.

   \return
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR 1 if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX ClearIVector (GRID *g)
{
  VECTOR *v;

  for (v = PFIRSTVECTOR(g); v != NULL; v = SUCCVC(v))
    VINDEX(v) = 0;

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Scale the interpolation matrix

   \param g - pointer to a grid
   \param theVD - vector descriptor

   This function scales all interpolation matrix entries by 1 / VINDEX
   resets the VINDEX for all vectors.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR 1 if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX ScaleIVector (GRID *g, VECDATA_DESC *theVD)
{
  VECTOR *v;
  INT i,comp;
  DOUBLE scale;

  i = 0;
  if (VD_IS_SCALAR(theVD))
  {
    comp = VD_SCALCMP(theVD);
    for (v = PFIRSTVECTOR(g); v != NULL; v = SUCCVC(v))
    {
      if (VINDEX(v) > 1)
        VVALUE(v,comp) *= (1.0 / VINDEX(v));
      VINDEX(v) = i++;
    }
    return (NUM_OK);
  }

  for (v = PFIRSTVECTOR(g); v != NULL; v = SUCCVC(v))
  {
    if (VINDEX(v) > 1)
    {
      scale = 1.0 / VINDEX(v);
      for (i=0; i<VD_NCMPS_IN_TYPE(theVD,VTYPE(v)); i++)
        VVALUE(v,VD_CMP_OF_TYPE(theVD,VTYPE(v),i)) *= scale;
    }
    VINDEX(v) = i++;
  }

  return (NUM_OK);
}



/****************************************************************************/
/** \brief Get the local interpolation matrix

   \param theElement - pointer to an element
   \param theFather - pointer to the father element
   \param me - number of nodal values of the element
   \param IntMat - local interpolation matrix
   \param theVD - vector descriptor

   This function gets the local interpolation matrix from the global
   interpolation matrix.

   \return
   .n    m    number of values in the father element
 */
/****************************************************************************/

INT NS_DIM_PREFIX GetInterpolationMatrix (ELEMENT *theElement, ELEMENT *theFather,
                                          INT me, DOUBLE *IntMat, VECDATA_DESC *theVD)
{
  VECTOR *fvec[MAX_NODAL_VECTORS];
  VECTOR *evec[MAX_NODAL_VECTORS];
  MATRIX *m;
  INT nev,nfv,ie,jf,ke,kf,nce,ncf;
  register SHORT i,j;
  DOUBLE *mptr;

  nev = GetAllVectorsOfElementOfType (theElement,evec,theVD);
  nfv = GetAllVectorsOfElementOfType (theFather,fvec,theVD);

  if (VD_IS_SCALAR(theVD)) {
    for (ie=0; ie<nev; ie++)
      for (jf=0; jf<nfv; jf++) {
        m = GetIMatrix(evec[ie],fvec[jf]);
        if (m == NULL)
          IntMat[jf*me+ie] = 0.0;
        else
          IntMat[jf*me+ie] = MVALUE(m,0);
      }
    return (nfv);
  }
  ke = 0;
  for (ie=0; ie<nev; ie++) {
    nce = VD_NCMPS_IN_TYPE(theVD,VTYPE(evec[ie]));
    kf = 0;
    for (jf=0; jf<nfv; jf++) {
      ncf = VD_NCMPS_IN_TYPE(theVD,VTYPE(fvec[jf]));
      m = GetIMatrix(evec[ie],fvec[jf]);
      if (m == NULL) {
        for (i=0; i<ncf; i++)
          for (j=0; j<nce; j++)
            IntMat[(kf+i)*me+ke+j] = 0.0;
      }
      else {
        mptr = MVALUEPTR(m,0);
        for (i=0; i<ncf; i++)
          for (j=0; j<nce; j++)
            IntMat[(kf+i)*me+ke+j] = mptr[i*nce+j];
      }
      kf += ncf;
    }
    ke += nce;
  }
  return(ncf);
}


static INT CheckBlock(INT me, INT ke, INT kf, INT ncf, INT nce, DOUBLE *IntMat)
{
  DOUBLE s;
  INT i,j;

  if (ncf != nce) {
    for (i=0; i<ncf; i++)
      for (j=0; j<nce; j++)
        if (IntMat[(kf+i)*me+ke+j] != 0.0)
          return(FULL_BLOCK);
    return(ZERO_BLOCK);
  }
  else {
    for (i=0; i<ncf; i++)
      for (j=0; j<nce; j++)
        if (i != j)
          if (IntMat[(kf+i)*me+ke+j] != 0.0)
            return(FULL_BLOCK);
    s = IntMat[kf*me+ke];
    if (s == 0.0) {
      for (i=1; i<ncf; i++)
        if (IntMat[(kf+i)*me+ke+i] != 0.0)
          return(FULL_BLOCK);
      return(ZERO_BLOCK);
    }
    for (i=1; i<ncf; i++)
      if (IntMat[(kf+i)*me+ke+i] != s)
        return(FULL_BLOCK);

    return(DIAG_BLOCK);
  }
}

/****************************************************************************/
/** \brief Add the local interpolation matrix

   \param theGrid - pointer to a grid
   \param theElement - pointer to an element
   \param theFather - pointer to the father element
   \param me - number of nodal values of the element
   \param IntMat - local interpolation matrix
   \param theVD - vector descriptor

   This function adds the local interpolation matrix to the global
   interpolation matrix.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR 1 if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX AddInterpolationMatrix (GRID *theGrid,
                                          ELEMENT *theElement, ELEMENT *theFather,
                                          INT me, DOUBLE *IntMat, VECDATA_DESC *theVD)
{
  VECTOR *fvec[MAX_NODAL_VECTORS];
  VECTOR *evec[MAX_NODAL_VECTORS];
  MATRIX *m;
  INT nev,nfv,ie,jf,ke,kf,nce,ncf,block;
  register SHORT i,j;
  DOUBLE *mptr,val;

  nev = GetAllVectorsOfElementOfType (theElement,evec,theVD);
  nfv = GetAllVectorsOfElementOfType (theFather,fvec,theVD);

  if (VD_IS_SCALAR(theVD))
  {
    for (ie=0; ie<nev; ie++)
    {
      for (jf=0; jf<nfv; jf++)
      {
        val = IntMat[jf*me+ie];
        if (val != 0.0)
        {
          m = GetIMatrix(evec[ie],fvec[jf]);
          if (m == NULL)
          {
            m = CreateIMatrix (theGrid,evec[ie],fvec[jf]);
            if (m == NULL)
              return(NUM_ERROR);
          }
          MVALUE(m,0) += val;
        }
      }
      VINDEX(evec[ie])++;
    }
    return (NUM_OK);
  }

  ke = 0;
  for (ie=0; ie<nev; ie++)
  {
    nce = VD_NCMPS_IN_TYPE(theVD,VTYPE(evec[ie]));
    kf = 0;
    for (jf=0; jf<nfv; jf++)
    {
      ncf = VD_NCMPS_IN_TYPE(theVD,VTYPE(fvec[jf]));
      block = CheckBlock(me,ke,kf,ncf,nce,IntMat);
      if (block != ZERO_BLOCK)
      {
        m = GetIMatrix(evec[ie],fvec[jf]);
        if (m == NULL)
        {
          m = CreateIMatrix (theGrid,evec[ie],fvec[jf]);
          if (m == NULL)
            return(NUM_ERROR);
        }
        if (block == DIAG_BLOCK) SETMDIAG(m,1);
        mptr = MVALUEPTR(m,0);
        for (i=0; i<ncf; i++)
          for (j=0; j<nce; j++)
            mptr[i*nce+j] += IntMat[(kf+i)*me+ke+j];
      }
      kf += ncf;
    }
    VINDEX(evec[ie])++;
    ke += nce;
  }

  return(NUM_OK);
}

/******** checks if damp != 1.0 ***********/

static INT CheckDamp (INT n, const DOUBLE *damp)
{
  INT i;

  for (i=0; i<n; i++)
    if (damp[i] != 1.0)
      return(1);

  return(0);
}


static INT RestrictByMatrix_General (GRID *FineGrid, const VECDATA_DESC *to,
                                     const VECDATA_DESC *from, const DOUBLE *damp, int flag)
{
  MATRIX *m;
  VECTOR *v,*w;
  DOUBLE sum,*vptr,*wptr,*mptr;
  INT vtype,wtype,vncomp,wncomp,vecskip;
  register SHORT i,j,xc,yc,xmask,ymask;
  const SHORT *offset;
  short rcomp;

  if (DOWNGRID(FineGrid)==NULL)
    return (NUM_NO_COARSER_GRID);

  if (VD_IS_SCALAR(to) && VD_IS_SCALAR(from))
  {
    xc    = VD_SCALCMP(to);
    yc    = VD_SCALCMP(from);
    xmask = VD_SCALTYPEMASK(to);
    ymask = VD_SCALTYPEMASK(from);

    for (w=PFIRSTVECTOR(DOWNGRID(FineGrid)); w!= NULL; w=SUCCVC(w))
      if ( (VDATATYPE(w)&xmask) && (VNCLASS(w)>=NEWDEF_CLASS) )
        VVALUE(w,xc) = 0.0;

    if (flag&1)
      rcomp = 1;
    else
      rcomp = 0;

    for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v))
      if ( (VDATATYPE(v)&ymask) && (VCLASS(v)>=NEWDEF_CLASS) )
        for (m=VISTART(v); m!= NULL; m = NEXT(m))
        {
          w = MDEST(m);
          if ( (VDATATYPE(w)&xmask) && (VECSKIP(w) == 0)
               && (CRITBIT(v, 0) == 0) )
            VVALUE(w,xc) += MVALUE(m,rcomp) * VVALUE(v,yc);
        }
    if (damp[0] != 1.0)
      for (w=PFIRSTVECTOR(DOWNGRID(FineGrid)); w!= NULL; w=SUCCVC(w))
        if ( (VDATATYPE(w)&xmask) && (VNCLASS(w)>=NEWDEF_CLASS) )
          VVALUE(w,xc) *= damp[0];

    return (NUM_OK);
  }

  for (w=PFIRSTVECTOR(DOWNGRID(FineGrid)); w!= NULL; w=SUCCVC(w))
    if (VNCLASS(w)>=NEWDEF_CLASS)
    {
      wtype = VTYPE(w);
      wncomp = VD_NCMPS_IN_TYPE(to,wtype);
      wptr = VVALUEPTR(w,VD_CMP_OF_TYPE(to,wtype,0));
      for (i=0; i<wncomp; i++)
        wptr[i] = 0.0;
    }

  for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v))
  {
    if (VCLASS(v)<NEWDEF_CLASS)
      continue;
    vtype = VTYPE(v);
    vncomp = VD_NCMPS_IN_TYPE(from,vtype);
    vptr = VVALUEPTR(v,VD_CMP_OF_TYPE(from,vtype,0));
    for (m=VISTART(v); m!= NULL; m = NEXT(m))
    {
      w = MDEST(m);
      mptr = MVALUEPTR(m,0);
      wtype = VTYPE(w);
      vecskip = VECSKIP(w);
      wncomp = VD_NCMPS_IN_TYPE(to,wtype);
      wptr = VVALUEPTR(w,VD_CMP_OF_TYPE(to,wtype,0));
      if (flag&1)
        rcomp=wncomp*vncomp;
      else
        rcomp=0;

      if (vecskip == 0)
        for (i=0; i<wncomp; i++)
        {
          sum = 0.0;
          for (j=0; j<vncomp; j++)
            sum += mptr[rcomp++] * vptr[j];
          wptr[i] += sum;
        }
      else
        for (i=0; i<wncomp; i++)
        {
          if (!(vecskip & (1<<i)))
          {
            sum = 0.0;
            for (j=0; j<vncomp; j++)
            {
              if(CRITBIT(v, j))
                rcomp++;
              else
                sum += mptr[rcomp++] * vptr[j];
            }
            wptr[i] += sum;
          }
          else
            rcomp += vncomp;
        }
    }
  }

  if (CheckDamp(VD_NCOMP(to),damp))
  {
    offset = VD_OFFSETPTR(to);
    for (w=PFIRSTVECTOR(DOWNGRID(FineGrid)); w!= NULL; w=SUCCVC(w))
      if (VNCLASS(w)>=NEWDEF_CLASS)
      {
        wtype = VTYPE(w);
        wncomp = VD_NCMPS_IN_TYPE(to,wtype);
        wptr = VVALUEPTR(w,VD_CMP_OF_TYPE(to,wtype,0));
        for (i=0; i<wncomp; i++)
          wptr[i] *= damp[offset[wtype]+i];
      }
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Restrict defect of fine vectors with NEWDEFECT_CLASS

   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor
   \param damp - damping factor for every component

   This function restricts defect of fine vectors with NEWDEFECT_CLASS,
   considers the VECSKIP-flags. It uses the transposed of the assembled
   interpolation matrix.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX RestrictByMatrix (GRID *FineGrid, const VECDATA_DESC *to,
                                    const VECDATA_DESC *from, const DOUBLE *damp)
{
  return( RestrictByMatrix_General (FineGrid,to,from,damp,0) );
}

INT NS_DIM_PREFIX RestrictByMatrix_s (GRID *FineGrid, const VECDATA_DESC *to,
                                      const VECDATA_DESC *from, const DOUBLE *damp)
{
  return( RestrictByMatrix_General (FineGrid,to,from,damp,1) );
}

/****************************************************************************/
/** \brief
   InterpolateCorrectionByMatrix - Interpolates the correction of the coarse grids

   SYNOPSIS:
   INT InterpolateCorrectionByMatrix (GRID *FineGrid, const VECDATA_DESC *to,
   const VECDATA_DESC *from, const DOUBLE *damp);

   PARAMETERS:
   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor
   \param damp - damping factor for every component

   DESCRIPTION:
   This function interpolates correction from coarse edge vectors,
   considers the VECSKIP-flags.
   It uses the transposed of the assembled interpolation matrix.

   \return
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

static INT InterpolateCorrectionByMatrix_General (GRID *FineGrid, const VECDATA_DESC *to,
                                                  const VECDATA_DESC *from,
                                                  const DOUBLE *damp, INT flag)
{
  MATRIX *m;
  VECTOR *v,*w;
  DOUBLE sum,*vptr,*wptr,*mptr;
  INT vtype,wtype,vncomp,wncomp,vecskip;
  register SHORT i,j,xc,yc,xmask,ymask;

  if (DOWNGRID(FineGrid)==NULL)
    return (NUM_NO_COARSER_GRID);

  dset(MYMG(FineGrid),GLEVEL(FineGrid),GLEVEL(FineGrid),ALL_VECTORS,to,0.0);

  if (VD_IS_SCALAR(to) && VD_IS_SCALAR(from))
  {
    xc    = VD_SCALCMP(to);
    yc    = VD_SCALCMP(from);
    xmask = VD_SCALTYPEMASK(to);
    ymask = VD_SCALTYPEMASK(from);

    for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v))
    {
      if ( (VDATATYPE(v)&xmask) )
      {
        if ((flag&1)==0)
          if (VECSKIPBIT(v, 0) ||
              CRITBIT(v, 0))
            continue;

        for (m=VISTART(v); m!= NULL; m = NEXT(m))
        {
          w = MDEST(m);
          if ( (VDATATYPE(w)&ymask) )
            VVALUE(v,xc) += MVALUE(m,0) * VVALUE(w,yc);
        }
      }
    }
    if (damp[0] != 1.0)
      if (dscalx(MYMG(FineGrid),GLEVEL(FineGrid),GLEVEL(FineGrid),
                 ALL_VECTORS,to,(DOUBLE*)damp))
        return (NUM_ERROR);

    return (NUM_OK);
  }

  for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v))
  {
    vtype = VTYPE(v);
    vncomp = VD_NCMPS_IN_TYPE(to,vtype);
    vptr = VVALUEPTR(v,VD_CMP_OF_TYPE(to,vtype,0));
    vecskip = VECSKIP(v);
    if (vecskip == 0 || ((flag&1)==1) )
    {
      for (m=VISTART(v); m!= NULL; m = NEXT(m))
      {
        w = MDEST(m);
        mptr = MVALUEPTR(m,0);
        wtype = VTYPE(w);
        wptr = VVALUEPTR(w,VD_CMP_OF_TYPE(from,wtype,0));
        wncomp = VD_NCMPS_IN_TYPE(from,wtype);
        for (i=0; i<vncomp; i++)
        {
          sum = 0.0;
          for (j=0; j<wncomp; j++)
            sum += mptr[j*vncomp+i] * wptr[j];
          vptr[i] += sum;
        }
      }
    }
    else
    {
      for (m=VISTART(v); m!= NULL; m = NEXT(m))
      {
        w = MDEST(m);
        mptr = MVALUEPTR(m,0);
        wtype = VTYPE(w);
        wptr = VVALUEPTR(w,VD_CMP_OF_TYPE(from,wtype,0));
        wncomp = VD_NCMPS_IN_TYPE(from,wtype);
        for (i=0; i<vncomp; i++)
          if (!(vecskip & (1<<i)) && !CRITBIT(v,i))
          {
            sum = 0.0;
            for (j=0; j<wncomp; j++)
              sum += mptr[j*vncomp+i] * wptr[j];
            vptr[i] += sum;
          }
      }
    }
  }

  if (CheckDamp(VD_NCOMP(to),damp))
    if (dscalx(MYMG(FineGrid),GLEVEL(FineGrid),GLEVEL(FineGrid),
               ALL_VECTORS,to,(DOUBLE*)damp))
      return (NUM_ERROR);

  return (NUM_OK);
}

INT NS_DIM_PREFIX InterpolateCorrectionByMatrix (GRID *FineGrid, const VECDATA_DESC *to,
                                                 const VECDATA_DESC *from,
                                                 const DOUBLE *damp)
{
  return (InterpolateCorrectionByMatrix_General(FineGrid,to,from,damp,0));
}

INT NS_DIM_PREFIX InterpolateCorrectionByMatrix_NoSkip (GRID *FineGrid, const VECDATA_DESC *to,
                                                        const VECDATA_DESC *from,
                                                        const DOUBLE *damp)
{
  return (InterpolateCorrectionByMatrix_General (FineGrid,to,from,damp,1));
}


/****************************************************************************/
/** \brief Interpolates the correction of the coarse grids

   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor
   \param damp - damping factor for every component

   This function interpolates correction from coarse edge vectors,
   considers the VECSKIP-flags.
   It uses the transposed of the assembled interpolation matrix.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX InterpolateNewVectorsByMatrix (GRID *FineGrid, const VECDATA_DESC *sol)
{
  MATRIX *m;
  VECTOR *v,*w;
  DOUBLE sum,*vptr,*wptr,*mptr;
  INT vtype,wtype,vncomp,wncomp,vecskip;
  register SHORT i,j,xc,xmask;

  if (DOWNGRID(FineGrid)==NULL)
    return (NUM_NO_COARSER_GRID);

  if (VD_IS_SCALAR(sol))
  {
    xc    = VD_SCALCMP(sol);
    xmask = VD_SCALTYPEMASK(sol);

    for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v))
      if ( (VDATATYPE(v)&xmask) && (VECSKIP(v) == 0) && (VNEW(v)) )
      {
        VVALUE(v,xc) = 0.0;
        for (m=VISTART(v); m!= NULL; m = NEXT(m))
        {
          w = MDEST(m);
          if ( (VDATATYPE(w)&xmask) )
            VVALUE(v,xc) += MVALUE(m,0) * VVALUE(w,xc);
        }
      }

    return (NUM_OK);
  }

  for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v))
  {
    if (!VNEW(v))
      continue;
    vtype = VTYPE(v);
    vncomp = VD_NCMPS_IN_TYPE(sol,vtype);
    vptr = VVALUEPTR(v,VD_CMP_OF_TYPE(sol,vtype,0));
    vecskip = VECSKIP(v);

    vecskip = 0;

    if (vecskip == 0)
    {
      for (i=0; i<vncomp; i++)
        vptr[i] = 0.0;
      for (m=VISTART(v); m!= NULL; m = NEXT(m))
      {
        w = MDEST(m);
        mptr = MVALUEPTR(m,0);
        wtype = VTYPE(w);
        wptr = VVALUEPTR(w,VD_CMP_OF_TYPE(sol,wtype,0));
        wncomp = VD_NCMPS_IN_TYPE(sol,wtype);
        for (i=0; i<vncomp; i++)
        {
          sum = 0.0;
          for (j=0; j<wncomp; j++)
            sum += mptr[j*vncomp+i] * wptr[j];
          vptr[i] += sum;
        }
      }
    }
    else
    {
      for (i=0; i<vncomp; i++)
        if (!(vecskip & (1<<i)))
          vptr[i] = 0.0;
      for (m=VISTART(v); m!= NULL; m = NEXT(m))
      {
        w = MDEST(m);
        mptr = MVALUEPTR(m,0);
        wtype = VTYPE(w);
        wptr = VVALUEPTR(w,VD_CMP_OF_TYPE(sol,wtype,0));
        wncomp = VD_NCMPS_IN_TYPE(sol,wtype);
        for (i=0; i<vncomp; i++)
          if (!(vecskip & (1<<i)))
          {
            sum = 0.0;
            for (j=0; j<wncomp; j++)
              sum += mptr[j*vncomp+i] * wptr[j];
            vptr[i] += sum;
          }
      }
    }
  }

  return (NUM_OK);
}


INT NS_DIM_PREFIX ClearGhostMatrix (GRID *g, MATDATA_DESC *Mat)
{
  VECTOR *v,*w;
  MATRIX *m;
  register DOUBLE *mptr;
  INT vtype,mtype,rmask,cmask,mc,vncomp,wncomp;
  register SHORT i,j,*mcomp;


  if (MD_IS_SCALAR(Mat)) {
    mc = MD_SCALCMP(Mat);
    rmask = MD_SCAL_RTYPEMASK(Mat);
    cmask = MD_SCAL_CTYPEMASK(Mat);

    for (v=PFIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
      if (VDATATYPE(v)&rmask)
        for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
          w = MDEST(m);
          if (!(VDATATYPE(w)&cmask))
            continue;
          MVALUE(m,mc) = 0.0;
        }
    return (NUM_OK);
  }
  for (v=PFIRSTVECTOR(g); v!=NULL; v=SUCCVC(v)) {
    vtype = VTYPE(v);
    for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
      w = MDEST(m);
      mtype = MTP(vtype,VTYPE(w));
      vncomp = MD_ROWS_IN_MTYPE(Mat,mtype);
      if (vncomp == 0) continue;
      wncomp = MD_COLS_IN_MTYPE(Mat,mtype);
      if (wncomp == 0) continue;
      mcomp = MD_MCMPPTR_OF_MTYPE(Mat,MTP(vtype,VTYPE(w)));
      mptr = MVALUEPTR(m,0);
      for (i=0; i<vncomp; i++)
        for (j=0; j<wncomp; j++)
          mptr[mcomp[i*wncomp+j]] = 0.0;
    }
  }

  return (NUM_OK);
}
/****************************************************************************/
/** \brief Galerkin assembling of the stiffness matrix

   \param FineGrid - pointer to grid
   \param Mat - matrix descriptor
   \param symmetric - flag for symmetry

   This function computes the Galerkin stiffness matrix for the given
   stiffness matrix on the fine grid and the interpolation matrix.
   The restriction matrix is the transposed matrix of the interpolation.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX AssembleGalerkinByMatrix (GRID *FineGrid, MATDATA_DESC *Mat, INT symmetric)
{
  GRID *CoarseGrid;
  MATRIX *m,*im,*jm,*cm;
  VECTOR *v,*w,*iv,*jv;
  register DOUBLE sum,*mptr,*cmptr,*imptr,*jmptr,mvalue;
  register DOUBLE *imptr0,*imptr1,*jmptr0,*jmptr1,*madjptr,fac;
  INT vtype,ivtype,mtype,cmtype,vncomp,wncomp,ivncomp,jvncomp,rmask,cmask;
  INT wtype,ivindex;
  register SHORT i,j,k,l,mc,*mcomp,*cmcomp,*madjcomp;

  CoarseGrid = DOWNGRID(FineGrid);
  if (CoarseGrid == NULL)
    return (NUM_NO_COARSER_GRID);

    #ifdef ModelP
  ClearGhostMatrix(CoarseGrid,Mat);
    #endif

  if (MD_IS_SCALAR(Mat)) {
    mc = MD_SCALCMP(Mat);
    rmask = MD_SCAL_RTYPEMASK(Mat);
    cmask = MD_SCAL_CTYPEMASK(Mat);

    for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v))
      if (VDATATYPE(v)&rmask)
        for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
          w = MDEST(m);
          if (!(VDATATYPE(w)&cmask))
            continue;
          mvalue = MVALUE(m,mc);
          for (im=VISTART(v); im!= NULL; im = NEXT(im)) {
            iv = MDEST(im);
            if (!(VDATATYPE(iv)&rmask))
              continue;
            fac = mvalue*MVALUE(im,0);
            for (jm=VISTART(w); jm!= NULL; jm = NEXT(jm)) {
              jv = MDEST(jm);
              if (!(VDATATYPE(jv)&cmask))
                continue;
              cm = GetMatrix(iv,jv);
              if (cm == NULL)
                cm = CreateExtraConnection(CoarseGrid,iv,jv);
              if (cm !=NULL)
                MVALUE(cm,mc) +=
                  fac * MVALUE(jm,0);
              else {                                           /* connection not in pattern */
                cm = GetMatrix(iv,iv);
                ASSERT(cm != NULL);
                MVALUE(cm,mc) +=
                  fac * MVALUE(jm,0);
              }
            }
          }
        }
  }
  else if (!symmetric) {
    for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
        w = MDEST(m);
        mtype = MTP(vtype,VTYPE(w));
        vncomp = MD_ROWS_IN_MTYPE(Mat,mtype);
        if (vncomp == 0) continue;
        wncomp = MD_COLS_IN_MTYPE(Mat,mtype);
        if (wncomp == 0) continue;
        mcomp = MD_MCMPPTR_OF_MTYPE(Mat,MTP(vtype,VTYPE(w)));
        mptr = MVALUEPTR(m,0);
        for (im=VISTART(v); im!= NULL; im = NEXT(im)) {
          iv = MDEST(im);
          ivtype = VTYPE(iv);
          imptr = MVALUEPTR(im,0);
          for (jm=VISTART(w); jm!= NULL; jm = NEXT(jm)) {
            jv = MDEST(jm);
            GET_MATRIX(iv,jv,cm);
            if (cm == NULL)
              cm = CreateExtraConnection(CoarseGrid,iv,jv);
            if (cm !=NULL) {
              cmtype = MTP(ivtype,VTYPE(jv));
              cmptr = MVALUEPTR(cm,0);
              cmcomp = MD_MCMPPTR_OF_MTYPE(Mat,cmtype);
              jmptr = MVALUEPTR(jm,0);
              ivncomp = MD_ROWS_IN_MTYPE(Mat,cmtype);
              jvncomp = MD_COLS_IN_MTYPE(Mat,cmtype);
              imptr1 = imptr;
              for (i=0; i<ivncomp; i++) {
                jmptr1 = jmptr;
                for (j=0; j<jvncomp; j++) {
                  sum = 0.0;
                  imptr0 = imptr1;
                  for (k=0; k<vncomp; k++) {
                    jmptr0 = jmptr1;
                    for (l=0; l<wncomp; l++) {
                      sum += *imptr0
                             * mptr[mcomp[k*wncomp+l]]
                             * *jmptr0++;
                    }
                    imptr0++;
                  }
                  cmptr[cmcomp[i*jvncomp+j]] += sum;
                  jmptr1 += wncomp;
                }
                imptr1 += vncomp;
              }
            }
            else {                                     /* connection not in pattern */
              cm = VSTART(jv);
              ASSERT(cm !=NULL);
              cmtype = MTP(VTYPE(jv),VTYPE(jv));
              jvncomp = MD_ROWS_IN_MTYPE(Mat,cmtype);
              cm = GetMatrix(iv,iv);
              ASSERT(cm !=NULL);
              cmtype = MTP(ivtype,ivtype);
              cmptr = MVALUEPTR(cm,0);
              cmcomp = MD_MCMPPTR_OF_MTYPE(Mat,cmtype);
              jmptr = MVALUEPTR(jm,0);
              ivncomp = MD_ROWS_IN_MTYPE(Mat,cmtype);
              for (i=0; i<ivncomp; i++)
                for (j=0; j<jvncomp; j++) {
                  sum = 0.0;
                  for (k=0; k<vncomp; k++)
                    for (l=0; l<wncomp; l++) {
                      sum += imptr[i*vncomp+k]
                             * mptr[mcomp[k*wncomp+l]]
                             * jmptr[j*wncomp+l];
                    }
                  cmptr[cmcomp[i*ivncomp+i]] += sum;
                }
            }
          }
        }
      }
    }
  }
  else if (MD_SUCC_COMP(Mat)) {
    for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
        w = MDEST(m);
        mtype = MTP(vtype,VTYPE(w));
        vncomp = MD_ROWS_IN_MTYPE(Mat,mtype);
        if (vncomp == 0) continue;
        wncomp = MD_COLS_IN_MTYPE(Mat,mtype);
        if (wncomp == 0) continue;
        mptr=MVALUEPTR(m,MD_MCMP_OF_MTYPE(Mat,MTP(vtype,VTYPE(w)),0));
        for (im=VISTART(v); im!= NULL; im = NEXT(im)) {
          iv = MDEST(im);
          ivindex = VINDEX(iv);
          ivtype = VTYPE(iv);
          imptr = MVALUEPTR(im,0);
          for (jm=VISTART(w); jm!= NULL; jm = NEXT(jm)) {
            jv = MDEST(jm);
            if (iv == jv) {
              cm = VSTART(iv);
              if (cm == NULL) {
                cm = CreateExtraConnection(CoarseGrid,iv,jv);
                if (cm == NULL)
                  REP_ERR_RETURN (NUM_ERROR);
              }
              cmtype = MTP(ivtype,ivtype);
              cmptr = MVALUEPTR(cm,
                                MD_MCMP_OF_MTYPE(Mat,cmtype,0));
              jmptr = MVALUEPTR(jm,0);
              if (MDIAG(im) && MDIAG(jm)) {
                for (k=0; k<vncomp; k++)
                  for (l=0; l<wncomp; l++) {
                    cmptr[k*wncomp+l] +=
                      *imptr * mptr[k*wncomp+l] * *jmptr;
                  }
              }
              else {
                ivncomp = MD_ROWS_IN_MTYPE(Mat,cmtype);
                imptr1 = imptr;
                for (i=0; i<ivncomp; i++) {
                  jmptr1 = jmptr;
                  for (j=0; j<ivncomp; j++) {
                    sum = 0.0;
                    imptr0 = imptr1;
                    for (k=0; k<vncomp; k++) {
                      jmptr0 = jmptr1;
                      for (l=0; l<wncomp; l++) {
                        sum += *imptr0
                               * mptr[k*wncomp+l]
                               * *jmptr0++;
                      }
                      imptr0++;
                    }
                    cmptr[i*ivncomp+j] += sum;
                    jmptr1 += wncomp;
                  }
                  imptr1 += vncomp;
                }
              }
            }
            else {
              if (ivindex < VINDEX(jv)) continue;
              for (cm = VSTART(iv); cm != NULL; cm = NEXT(cm))
                if (MDEST(cm) == jv) break;
              if (cm == NULL) {
                cm = CreateExtraConnection(CoarseGrid,iv,jv);
                if (cm == NULL)
                  REP_ERR_RETURN (NUM_ERROR);
              }
              cmtype = MTP(ivtype,VTYPE(jv));
              cmptr = MVALUEPTR(cm,
                                MD_MCMP_OF_MTYPE(Mat,cmtype,0));
              jmptr = MVALUEPTR(jm,0);
              if (MDIAG(im) && MDIAG(jm)) {
                for (k=0; k<vncomp; k++)
                  for (l=0; l<wncomp; l++) {
                    cmptr[k*wncomp+l] +=
                      *imptr * mptr[k*wncomp+l] * *jmptr;
                  }
              }
              else {
                ivncomp = MD_ROWS_IN_MTYPE(Mat,cmtype);
                jvncomp = MD_COLS_IN_MTYPE(Mat,cmtype);
                imptr1 = imptr;
                for (i=0; i<ivncomp; i++) {
                  jmptr1 = jmptr;
                  for (j=0; j<jvncomp; j++) {
                    sum = 0.0;
                    imptr0 = imptr1;
                    for (k=0; k<vncomp; k++) {
                      jmptr0 = jmptr1;
                      for (l=0; l<wncomp; l++) {
                        sum += *imptr0
                               * mptr[k*wncomp+l]
                               * *jmptr0++;
                      }
                      imptr0++;
                    }
                    cmptr[i*jvncomp+j] += sum;
                    jmptr1 += wncomp;
                  }
                  imptr1 += vncomp;
                }
              }
            }
          }
        }
      }
    }
    for (v=PFIRSTVECTOR(CoarseGrid); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      ivindex = VINDEX(v);
      for (cm=VSTART(v); cm!=NULL; cm=MNEXT(cm)) {
        w = MDEST(cm);
        if (ivindex < VINDEX(w)) continue;
        wtype = MDESTTYPE(cm);
        mtype = MTP(vtype,wtype);
        mptr = MVALUEPTR(cm,MD_MCMP_OF_MTYPE(Mat,mtype,0));
        ivncomp = MD_ROWS_IN_MTYPE(Mat,mtype);
        jvncomp = MD_COLS_IN_MTYPE(Mat,mtype);
        madjptr = MVALUEPTR(MADJ(cm),
                            MD_MCMP_OF_MTYPE(Mat,MTP(wtype,vtype),0));
        for (i=0; i<ivncomp; i++)
          for (j=0; j<jvncomp; j++)
            madjptr[j*ivncomp+i] =
              mptr[i*jvncomp+j];
      }
    }
  }
  else if (symmetric) {
    for (v=FIRSTVECTOR(FineGrid); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
        w = MDEST(m);
        mtype = MTP(vtype,VTYPE(w));
        vncomp = MD_ROWS_IN_MTYPE(Mat,mtype);
        if (vncomp == 0) continue;
        wncomp = MD_COLS_IN_MTYPE(Mat,mtype);
        if (wncomp == 0) continue;
        mcomp = MD_MCMPPTR_OF_MTYPE(Mat,MTP(vtype,VTYPE(w)));
        mptr = MVALUEPTR(m,0);
        for (im=VISTART(v); im!= NULL; im = NEXT(im)) {
          iv = MDEST(im);
          ivindex = VINDEX(iv);
          ivtype = VTYPE(iv);
          imptr = MVALUEPTR(im,0);
          for (jm=VISTART(w); jm!= NULL; jm = NEXT(jm)) {
            jv = MDEST(jm);
            if (ivindex < VINDEX(jv)) continue;
            GET_MATRIX(iv,jv,cm);
            if (cm == NULL)
              cm = CreateExtraConnection(CoarseGrid,iv,jv);
            if (cm !=NULL) {
              cmtype = MTP(ivtype,VTYPE(jv));
              cmptr = MVALUEPTR(cm,0);
              cmcomp = MD_MCMPPTR_OF_MTYPE(Mat,cmtype);
              jmptr = MVALUEPTR(jm,0);
              ivncomp = MD_ROWS_IN_MTYPE(Mat,cmtype);
              jvncomp = MD_COLS_IN_MTYPE(Mat,cmtype);
              imptr1 = imptr;
              for (i=0; i<ivncomp; i++) {
                jmptr1 = jmptr;
                for (j=0; j<jvncomp; j++) {
                  sum = 0.0;
                  imptr0 = imptr1;
                  for (k=0; k<vncomp; k++) {
                    jmptr0 = jmptr1;
                    for (l=0; l<wncomp; l++) {
                      sum += *imptr0
                             * mptr[mcomp[k*wncomp+l]]
                             * *jmptr0++;
                    }
                    imptr0++;
                  }
                  cmptr[cmcomp[i*jvncomp+j]] += sum;
                  jmptr1 += wncomp;
                }
                imptr1 += vncomp;
              }
            }
            else {                                     /* connection not in pattern */
              cm = VSTART(jv);
              ASSERT(cm !=NULL);
              cmtype = MTP(VTYPE(jv),VTYPE(jv));
              jvncomp = MD_ROWS_IN_MTYPE(Mat,cmtype);
              cm = GetMatrix(iv,iv);
              ASSERT(cm !=NULL);
              cmtype = MTP(ivtype,ivtype);
              cmptr = MVALUEPTR(cm,0);
              cmcomp = MD_MCMPPTR_OF_MTYPE(Mat,cmtype);
              jmptr = MVALUEPTR(jm,0);
              ivncomp = MD_ROWS_IN_MTYPE(Mat,cmtype);
              for (i=0; i<ivncomp; i++)
                for (j=0; j<jvncomp; j++) {
                  sum = 0.0;
                  for (k=0; k<vncomp; k++)
                    for (l=0; l<wncomp; l++) {
                      sum += imptr[i*vncomp+k]
                             * mptr[mcomp[k*wncomp+l]]
                             * jmptr[j*wncomp+l];
                    }
                  cmptr[cmcomp[i*ivncomp+i]] += sum;
                }
            }
          }
        }
      }
    }
    for (v=PFIRSTVECTOR(CoarseGrid); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      ivindex = VINDEX(v);
      for (cm=VSTART(v); cm!=NULL; cm=MNEXT(cm)) {
        if (ivindex <= VINDEX(MDEST(cm))) continue;
        wtype = MDESTTYPE(cm);
        mtype = MTP(vtype,wtype);
        mptr = MVALUEPTR(cm,0);
        mcomp = MD_MCMPPTR_OF_MTYPE(Mat,mtype);
        ivncomp = MD_ROWS_IN_MTYPE(Mat,mtype);
        jvncomp = MD_COLS_IN_MTYPE(Mat,mtype);
        madjptr = MVALUEPTR(MADJ(cm),0);
        madjcomp = MD_MCMPPTR_OF_MTYPE(Mat,MTP(wtype,vtype));
        for (i=0; i<ivncomp; i++)
          for (j=0; j<jvncomp; j++)
            madjptr[madjcomp[j*ivncomp+i]] =
              mptr[mcomp[i*jvncomp+j]];
      }
    }
  }
        #ifdef ModelP
  if (l_ghostmatrix_collect(CoarseGrid,Mat))
    return(NUM_ERROR);
        #endif

  return(NUM_OK);
}

/****************************************************************************/
/** \brief Restriction for diagonally scaled mg

   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor
   \param Amat - matrix to compute weights
   \param damp - damping factor for every component

   This function restricts defect of fine node vectors with NEWDEFECT_CLASS
   to the next coarser grid.
   First, it resets all components to zero. It considers the VECSKIP-flags.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

static INT ScaledMGRestrictNodeVector (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp)
{
  GRID *CoarseGrid;
  NODE *theNode;
  VECTOR *v,*vc;
  const SHORT *toComp,*fromComp;
  INT i,k,ncomp,vecskip;
  INT dt;
  MATRIX *im;

  CoarseGrid = DOWNGRID(FineGrid);

  toComp   = VD_ncmp_cmpptr_of_otype(to,NODEVEC,&ncomp);
  if (ncomp <= 0) return(NUM_ERROR);
  if (ncomp>MAX_SINGLE_VEC_COMP) return (NUM_BLOCK_TOO_LARGE);
  fromComp = VD_cmpptr_of_otype(from,NODEVEC);

  /* reset coarser defect at positions where a new defect is restricted */
  dt = VD_DATA_TYPES(to);
  for (v=PFIRSTVECTOR(CoarseGrid); v!= NULL; v=SUCCVC(v))
    if ((VOTYPE(v)==NODEVEC) && V_IN_DATATYPE(v,dt))
      if (VNCLASS(v)>=NEWDEF_CLASS)
        for (i=0; i<ncomp; i++)
          VVALUE(v,toComp[i]) = 0.0;

  /* compute contributions to all coarse node vectors */
  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    v = NVECTOR(theNode);
    if (VCLASS(v)<NEWDEF_CLASS) continue;
    if (!V_IN_DATATYPE(v,dt)) continue;

    for (im=VISTART(v); im!=NULL; im=MNEXT(im))
    {
      vc = MDEST(im);
      vecskip = VECSKIP(vc);
      for (i=0; i<ncomp; i++)
        if (!(vecskip & (1<<i)))
          for (k=0; k<ncomp; k++) VVALUE(vc,toComp[i]) += MVALUE(im,i*ncomp+k)*VVALUE(v,fromComp[k]);
    }
  }

  return (NUM_OK);
}


/****************************************************************************/
/** \brief Matrix dependent restriction of fine vectors with NEWDEFECT_CLASS

   \param FineGrid - pointer to grid
   \param to - type vector descriptor
   \param from  - type vector descriptor
   \param Mat - fine grid matrix
   \param damp - damping factor for every component

   This function restricts defect of fine vectors with NEWDEFECT_CLASS,
   considers the VECSKIP-flags.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX ScaledMGRestrict (GRID *FineGrid, const VECDATA_DESC *to, const VECDATA_DESC *from, const DOUBLE *damp)
{
  INT vtype,rv;
  const SHORT *offset;

  if (DOWNGRID(FineGrid)==NULL)
    return (NUM_NO_COARSER_GRID);

  offset = VD_OFFSETPTR(to);

  for (vtype=0; vtype<NVECTYPES; vtype++)
    if (VD_ISDEF_IN_TYPE(to,vtype))
      switch (GetUniqueOTypeOfVType(MGFORMAT(MYMG(FineGrid)),vtype))
      {
      case ELEMVEC :
        PrintErrorMessage('E',"MatDepRestrict","only node vector is implemented");
        return(NUM_ERROR);
      case NODEVEC :
        if ((rv=ScaledMGRestrictNodeVector(FineGrid,to,from,damp+offset[vtype]))!=NUM_OK)
          return (rv);
        break;
      case EDGEVEC :
        PrintErrorMessage('E',"MatDepRestrict","only node vector is implemented");
        return(NUM_ERROR);
      case SIDEVEC :
        PrintErrorMessage('E',"MatDepRestrict","only node vector is implemented");
        return(NUM_ERROR);
      default :
        UserWrite("not unique");
        return (NUM_ERROR);
      }

  return (NUM_OK);
}




/****************************************************************************/
/** \brief Compute restriction matrix for scaled mg

   \param FineGrid - pointer to fine grid equations
   \param Mat - matrix to be computed

   This function computes the modified restriction matrix used in diagonally
   scaled multigrid algorithm.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

#undef _LOCAL_DEBUG_

INT NS_DIM_PREFIX InstallScaledRestrictionMatrix (GRID *FineGrid, const MATDATA_DESC *Mat, DOUBLE cut)
{
  NODE *theNode;
  VECTOR *vf,*vc;
  INT i,j,k,n,l,ncomp,nc,vecskip,A,rdt;
  SHORT *comps;

  DOUBLE Dcoarseinv[MAX_SINGLE_VEC_COMP*MAX_SINGLE_VEC_COMP];
  DOUBLE Q[MAX_SINGLE_VEC_COMP*MAX_SINGLE_VEC_COMP];
  DOUBLE F[MAX_SINGLE_VEC_COMP*MAX_SINGLE_VEC_COMP];

  DOUBLE *Dfine,*Dcoarse;
  MATRIX *im;
  ELEMENT *theElement;
  VERTEX *theVertex;
  DOUBLE c[MAX_CORNERS_OF_ELEM],s;
        #ifdef _LOCAL_DEBUG_
  VECTOR *v;
  char buffer[128];
        #endif

  /* we handle only node vectors here ! */
  comps = MD_nr_nc_mcmpptr_of_ro_co(Mat,NODEVEC,NODEVEC,&ncomp,&nc);
  if (ncomp <= 0) REP_ERR_RETURN(1);
  if (ncomp>MAX_SINGLE_VEC_COMP) REP_ERR_RETURN(NUM_ERROR);

  /* check matrix format and get components */
  if (nc!=ncomp) REP_ERR_RETURN(1);

  A = comps[0];
  for (i=0; i<ncomp*ncomp; i++)
    if (comps[i]!=A+i)
    {
      PrintErrorMessage('E',"InstallRestrictionMatrix","matrix format incorrect");
      REP_ERR_RETURN(1);
    }

  /* compute contributions of fine node to coarse nodes */
  rdt = MD_ROW_DATA_TYPES(Mat);
  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    vf = NVECTOR(theNode);
    if (VCLASS(vf)<NEWDEF_CLASS) continue;
    if (!V_IN_DATATYPE(vf,rdt)) continue;

    /* fine grid diagonal block */
    Dfine = &(MVALUE(VSTART(vf),A));

    if (CORNERTYPE(theNode))             /* This node is also in the coarse grid */
    {
      vc = NVECTOR((NODE *)NFATHER(theNode));

      /* compute Q = D(vc)(D(vf))^{-1} */
      Dcoarse = &(MVALUE(VSTART(vc),0));
      if (InvertSmallBlock(ncomp,comps,Dcoarse,Dcoarseinv)!=NUM_OK)
      {
        UserWriteF("ncomp=%d, comps[0]=%d, Dcoarse=%f\n",ncomp,comps[0],*Dcoarse);
        REP_ERR_RETURN(1);
      }
      for (i=0; i<ncomp; i++)
        for (j=0; j<ncomp; j++) {
          Q[i*ncomp+j] = 0.0;
          for (k=0; k<ncomp; k++)
            Q[i*ncomp+j] += Dcoarseinv[i*ncomp+k]*Dfine[k*ncomp+j];
        }

      /* apply filter ! */
      for (i=0; i<ncomp; i++)
        for (j=0; j<ncomp; j++)
        {
          Q[i*ncomp+j] = MAX(0.0,MIN(cut,Q[i*ncomp+j]));
        }

      /* allocate restriction matrix entry */
      im = GetIMatrix(vf,vc);
      if (im==NULL) {
        im = CreateIMatrix(FineGrid,vf,vc);
        if (im==NULL) {
          UserWrite("Could not create interpolation matrix\n");
          REP_ERR_RETURN(1);
        }
      }
      for (i=0; i<ncomp*ncomp; i++) MVALUE(im,i) = Q[i];

    }
    else             /* This node is only in fine grid */
    {
      theVertex  = MYVERTEX(theNode);
      theElement = VFATHER(theVertex);
      n = CORNERS_OF_ELEM(theElement);
      GNs(n,LCVECT(theVertex),c);
      for (l=0; l<n; l++)
      {
        /* for all corners */
        vc = NVECTOR(CORNER(theElement,l));
        vecskip = VECSKIP(vc);

        /* copy fine and eliminate dirichlet nodes in fine matrix ! */
        for (i=0; i<ncomp*ncomp; i++) F[i] = Dfine[i];
        for (i=0; i<ncomp; i++)
          if ((vecskip & (1<<i)))
            for (j=0; j<ncomp; j++)
              if (i==j) F[i] = 1.0;else F[i] = 0.0;

        /* compute Q = D(vc)^{-1}(D(vf)) */
        Dcoarse = &(MVALUE(VSTART(vc),0));
        if (InvertSmallBlock(ncomp,comps,Dcoarse,Dcoarseinv)!=NUM_OK)
          return (__LINE__);
        for (i=0; i<ncomp; i++) {
          if (!(vecskip & (1<<i))) s=1;else s=0;
          for (j=0; j<ncomp; j++) {
            Q[i*ncomp+j] = 0.0;
            for (k=0; k<ncomp; k++)
              Q[i*ncomp+j] += s*Dcoarseinv[i*ncomp+k]*F[k*ncomp+j];
          }
        }

        /* apply filter ! */
        for (i=0; i<ncomp; i++)
          for (j=0; j<ncomp; j++)
          {
            Q[i*ncomp+j] = MAX(0.0,MIN(cut,Q[i*ncomp+j]));
          }

        /* allocate restriction matrix entry */
        im = GetIMatrix(vf,vc);
        if (im==NULL) {
          im = CreateIMatrix(FineGrid,vf,vc);
          if (im==NULL) {
            UserWrite("Could not create interpolation matrix\n");
            REP_ERR_RETURN(1);
          }
        }
        for (i=0; i<ncomp*ncomp; i++) MVALUE(im,i) = Q[i]*c[l];
      }
    }
  }

        #ifdef _LOCAL_DEBUG_
  UserWriteF("---- LEVEL %d ----\n",GLEVEL(FineGrid));
  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    UserWriteF("*NODE (%12.4g,%12.4g) %ld\n",(DOUBLE)XC(MYVERTEX(theNode)),
               (DOUBLE)YC(MYVERTEX(theNode)),ID(theNode));

    v = NVECTOR(theNode);
    for (im=VISTART(v); im!=NULL; im=MNEXT(im))
    {
      vc = MDEST(im);
      vecskip = VECSKIP(vc);
      for (i=0; i<ncomp; i++)
      {
        UserWriteF(" DEST (%12.4g,%12.4g) %5ld, ",
                   (DOUBLE)XC(MYVERTEX(VMYNODE(vc))),(DOUBLE)YC(MYVERTEX(VMYNODE(vc))),ID(VMYNODE(vc)));
        for (k=0; k<ncomp*ncomp; k++)
          UserWriteF(" %12.4E",MVALUE(im,k));
        UserWrite("\n");
      }
    }
  }
        #endif

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Scale system of equations by point-block-diagonal

   \param FineGrid - pointer to grid
   \param Mat - matrix
   \param rhs - right hand side

   Scales Ax=b to DAx=Db, where D is the inverse of the diagonal blocks of A.

   \return
   .n    NUM_OK if ok
   .n    NUM_ERROR if error occured.
 */
/****************************************************************************/

INT NS_DIM_PREFIX DiagonalScaleSystem (GRID *FineGrid, const MATDATA_DESC *Mat, const MATDATA_DESC *ConsMat, const VECDATA_DESC *rhs)
{
  NODE *theNode;
  INT A,ConsA,b,n,i,j,k;
  VECTOR *vi;
  MATRIX *mij;
  SHORT *comps,*p;
  SHORT *ConsComps;
  DOUBLE Dfineinv[MAX_SINGLE_VEC_COMP*MAX_SINGLE_VEC_COMP];
  DOUBLE Q[MAX_SINGLE_VEC_COMP*MAX_SINGLE_VEC_COMP];
  DOUBLE r[MAX_SINGLE_VEC_COMP];
  DOUBLE *Dfine,*bfine;
  INT ncomp,nc;

  /* check assumptions and get components */

  /* rhs */
  p    = VD_ncmp_cmpptr_of_otype(rhs,NODEVEC,&n);
  if (n <= 0) REP_ERR_RETURN(1);
  if (n>MAX_SINGLE_VEC_COMP) REP_ERR_RETURN (1);

  /* Mat */
  comps = MD_nr_nc_mcmpptr_of_ro_co(Mat,NODEVEC,NODEVEC,&ncomp,&nc);
  if (ncomp != n) REP_ERR_RETURN(1);
  if (nc    != n) REP_ERR_RETURN(1);

  /* ConsMat */
  ConsComps = MD_nr_nc_mcmpptr_of_ro_co(ConsMat,NODEVEC,NODEVEC,&ncomp,&nc);
  if (ncomp != n) REP_ERR_RETURN(1);
  if (nc    != n) REP_ERR_RETURN(1);

  /* check subsequence of components */
  b = p[0];
  for (i=0; i<n; i++)
    if (p[i]!=b+i)
    {
      PrintErrorMessage('E',"ScaleSystem","vector format incorrect");
      REP_ERR_RETURN(NUM_ERROR);
    }
  A = comps[0];
  for (i=0; i<n*n; i++)
    if (comps[i]!=A+i)
    {
      PrintErrorMessage('E',"ScaleSystem","matrix format incorrect");
      REP_ERR_RETURN(NUM_ERROR);
    }
  ConsA = ConsComps[0];
  for (i=0; i<n*n; i++)
    if (ConsComps[i]!=ConsA+i)
    {
      PrintErrorMessage('E',"ScaleSystem","cons matrix format incorrect");
      REP_ERR_RETURN(NUM_ERROR);
    }

  /* scale system by (consistent !) point block diagonal */
  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    /* get vector */
    vi = NVECTOR(theNode);

    /* invert consistent diagonal block */
    Dfine = &(MVALUE(VSTART(vi),0));
    if (InvertSmallBlock(n,ConsComps,Dfine,Dfineinv)!=NUM_OK) REP_ERR_RETURN (NUM_ERROR);

    /* multiply row from left */
    for (mij=VSTART(vi); mij!=NULL; mij=MNEXT(mij))
    {
      /*
            if (CEXTRA(mij)) continue;
       */
      Dfine = &(MVALUE(mij,A));
      for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
          Q[i*n+j] = 0.0;
          for (k=0; k<n; k++)
            Q[i*n+j] += Dfineinv[i*n+k]*Dfine[k*n+j];
        }
      for (i=0; i<n*n; i++) Dfine[i] = Q[i];
    }

    /* and the right hand side */
    bfine = &(VVALUE(vi,b));
    for (i=0; i<n; i++) {
      r[i] = 0;
      for (j=0; j<n; j++)
        r[i] += Dfineinv[i*n+j]*bfine[j];
    }
    for (i=0; i<n; i++) bfine[i] = r[i];

  }

  return(NUM_OK);
}


INT NS_DIM_PREFIX CreateStandardNodeRestProl (GRID *FineGrid, INT ncomp)
{
  NODE *theNode;
  VECTOR *vf,*vc;
  INT i,j,n,l;
  MATRIX *im;
  ELEMENT *theElement;
  VERTEX *theVertex;
  DOUBLE c[MAX_CORNERS_OF_ELEM];

  for (theNode=FIRSTNODE(FineGrid); theNode!= NULL; theNode=SUCCN(theNode))
  {
    vf = NVECTOR(theNode);
    if (CORNERTYPE(theNode))             /* This node is also in the coarse grid */
    {
      vc = NVECTOR((NODE *)NFATHER(theNode));
      /* allocate restriction matrix entry */
      im = GetIMatrix(vf,vc);
      if (im==NULL) {
        im = CreateIMatrix(FineGrid,vf,vc);
        if (im==NULL) {
          UserWrite("Could not create interpolation matrix\n");
          return(__LINE__);
        }
      }
      for (i=0; i<ncomp; i++)
        for (j=0; j<ncomp; j++)
          if (i == j) MVALUE(im,i*ncomp+j) = 1.0;
          else MVALUE(im,i*ncomp+j) = 0.0;
    }
    else             /* This node is only in fine grid */
    {
      theVertex  = MYVERTEX(theNode);
      theElement = VFATHER(theVertex);
      n = CORNERS_OF_ELEM(theElement);
      GNs(n,LCVECT(theVertex),c);
      for (l=0; l<n; l++)
      {
        if (c[l] == 0.0) continue;
        vc = NVECTOR(CORNER(theElement,l));
        im = GetIMatrix(vf,vc);
        if (im==NULL) {
          im = CreateIMatrix(FineGrid,vf,vc);
          if (im==NULL) {
            UserWrite("Could not create interpolation matrix\n");
            return(__LINE__);
          }
        }
        for (i=0; i<ncomp; i++)
          for (j=0; j<ncomp; j++)
            if (i == j) MVALUE(im,i*ncomp+j) = c[l];
            else MVALUE(im,i*ncomp+j) = 0.0;
      }
    }
  }
  return (NUM_OK);
}
#endif /* interpolation matrix is defined */
