// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      fegeom.c                                                          */
/*                                                                          */
/* Purpose:   geometry related data for a general element in the cfe scheme */
/*                                                                          */
/* Author:	  Peter Bastian                                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  fon: 0049-(0)711-685-7003										*/
/*			  fax: 0049-(0)711-685-7000										*/
/*																			*/
/* History:   02.07.96 begin, ug version 3.0								*/
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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "devices.h"
#include "enrol.h"
#include "compiler.h"
#include "misc.h"
#include "gm.h"
#include "ugenv.h"
#include "ugm.h"
#include "algebra.h"
#include "cmdint.h"
#include "commands.h"
#include "helpmsg.h"
#include "shapes.h"
#include "quadrature.h"

#include "fegeom.h"

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
#define CP(p,p1)            (p)[0] = (p1)[0]; (p)[1] = (p1)[1]
#define AVG2(p,p1,p2)         p[0] = 0.5*(p1[0]+p2[0]); p[1] = 0.5*(p1[1]+p2[1])
#define AVG3(p,p1,p2,p3)      p[0] = (p1[0]+p2[0]+p3[0])/3.0; p[1] = (p1[1]+p2[1]+p3[1])/3.0
#define AVG4(p,p1,p2,p3,p4)   p[0] = 0.25*(p1[0]+p2[0]+p3[0]+p4[0]); p[1] = 0.25*(p1[1]+p2[1]+p3[1]+p4[1])
#endif


#ifdef __THREEDIM__
#define CP(p,p1) p[0] = p1[0]; p[1] = p1[1],  p[2] = p1[2]
#define AVG2(p,p1,p2) p[0] = 0.5*(p1[0]+p2[0]); p[1] = 0.5*(p1[1]+p2[1]); p[2] = 0.5*(p1[2]+p2[2])
#define AVG3(p,p1,p2,p3) p[0] = (p1[0]+p2[0]+p3[0])/3.0; p[1] = (p1[1]+p2[1]+p3[1])/3.0;  p[2] = (p1[2]+p2[2]+p3[2])/3.0
#define AVG4(p,p1,p2,p3,p4) p[0] = 0.25*(p1[0]+p2[0]+p3[0]+p4[0]); p[1] = 0.25*(p1[1]+p2[1]+p3[1]+p4[1]); p[2] = 0.25*(p1[2]+p2[2]+p3[2]+p4[2])
#endif

#define Xi  ((DOUBLE)ip_local[0])
#define Eta ((DOUBLE)ip_local[1])
#define Mu  ((DOUBLE)ip_local[2])

#define J11 J[0][0]
#define J12 J[0][1]
#define J13 J[0][2]
#define J21 J[1][0]
#define J22 J[1][1]
#define J23 J[1][2]
#define J31 J[2][0]
#define J32 J[2][1]
#define J33 J[2][2]

#define X0  co_global[0][0]
#define X1  co_global[1][0]
#define X2  co_global[2][0]
#define X3  co_global[3][0]
#define X4  co_global[4][0]
#define X5  co_global[5][0]
#define X6  co_global[6][0]
#define X7  co_global[7][0]

#define Y0  co_global[0][1]
#define Y1  co_global[1][1]
#define Y2  co_global[2][1]
#define Y3  co_global[3][1]
#define Y4  co_global[4][1]
#define Y5  co_global[5][1]
#define Y6  co_global[6][1]
#define Y7  co_global[7][1]

#define Z0  co_global[0][2]
#define Z1  co_global[1][2]
#define Z2  co_global[2][2]
#define Z3  co_global[3][2]
#define Z4  co_global[4][2]
#define Z5  co_global[5][2]
#define Z6  co_global[6][2]
#define Z7  co_global[7][2]

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

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* definition of functions													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   EvaluateFEGeometry - compute geometry information for given element

   SYNOPSIS:
   INT EvaluateFEGeometry (ELEMENT *e, FEElementGeometry *eg);

   PARAMETERS:
   .  e - given element
   .  eg - data structure to be filled

   DESCRIPTION:
   This routine fills the given data structure with geometry related
   values

   RETURN VALUES:
   0 when o.k.

   1 if an error occured.
   D*/
/****************************************************************************/

INT EvaluateFEGeometry (ELEMENT *e, FEElementGeometry *geo)
{
  VERTEX *v;
  QUADRATURE *q;
  GaussPoint *gp;
  Edge *ed;
  BoundarySide *bs;
  BoundaryGaussPoint *bgp;
  int i,j,k;
  DOUBLE nvals[MAX_CORNERS_OF_ELEM];
  DOUBLE sco_global[MAX_CORNERS_OF_ELEM][DIM];
  DOUBLE refElemVolume;
  DOUBLE refFaceArea;

  /* general info */
  geo->e = e;
  geo->tag = TAG(e);
  geo->nc = CORNERS_OF_ELEM(e);
  geo->ned = 0;       /* initially */
  geo->nbs = 0;       /* initially */
  geo->ngp = 0;       /* initially */

  /* corners */
  for (i=0; i<geo->nc; i++) {
    v = MYVERTEX(CORNER(e,i));
    CP(geo->co_global[i],CVECT(v));
    LocalCornerCoordinates(DIM,geo->tag,i,geo->co_local[i]);
    geo->node_property[i] =  NPROP(CORNER(e,i));
  }

  /* the gauss points */
  if (DIM==2)       /* select a quadrature rule */
    switch (geo->nc)
    {
    case 3 :
      q = GetQuadrature(2,3,2);
      refElemVolume = 0.5;
      break;
    case 4 :
      q = GetQuadrature(2,4,4);
      refElemVolume = 1.0;
      break;
    default :
      return(1);
    }
  if (DIM==3)
    switch (geo->nc)
    {
    case 4 :
      q = GetQuadrature(3,4,2);
      refElemVolume = 1.0/6.0;
      break;
    case 5 :
      q = GetQuadrature(3,5,2);
      refElemVolume = 1.0/3.0;
      break;
    case 6 :
      q = GetQuadrature(3,6,2);
      refElemVolume = 1.0/2.0;
      break;
    case 8 :
      q = GetQuadrature(3,8,2);
      refElemVolume = 1.0;
      break;
    default :
      return(1);
    }
  if (q==NULL) return(1);
  geo->ngp = Q_NIP(q);
  for (i=0; i<geo->ngp; i++)
  {
    /* get gauss point */
    gp = geo->gp+i;

    /* get weight */
    gp->weight = refElemVolume*Q_WEIGHT(q,i);

    /* get local coordinates */
    for (j=0; j<DIM; j++) gp->local[j] = Q_LOCAL(q,i)[j];

    /* evaluate basis functions and gradient */
    for (j=0; j<geo->nc; j++)
    {
      gp->N[j] = GN(geo->nc,j,Q_LOCAL(q,i));
      D_GN(geo->nc,j,Q_LOCAL(q,i),gp->gradN[j]);
    }

    /* inverse of jacobian */
    JacobianInverse(DIM,geo->tag,geo->co_global,gp->local,gp->Jinv,&gp->AbsdetJ);
    gp->AbsdetJ = ABS(gp->AbsdetJ);

  }

  /* loop over all connections */
  for (i=0; i<geo->nc; i++)
    for (j=i+1; j<geo->nc; j++)
    {
      /* alloc a new connection */
      ed = geo->ed+geo->ned;
      geo->ned++;

      /* end points of connection */
      ed->i = i; ed->j = j;

      /* edge vector pointing from i to j */
      for (k=0; k<DIM; k++) ed->edge[k] = geo->co_global[j][k] - geo->co_global[i][k];

      /* edge midpoint in local coordinates */
      AVG2(ed->local,geo->co_local[i],geo->co_local[j]);

      /* inverse of jacobian */
      JacobianInverse(DIM,geo->tag,geo->co_global,ed->local,ed->Jinv,&ed->detJ);
      ed->detJ = ABS(ed->detJ);
    }

  /* the boundary */
  if (OBJT(e)==BEOBJ)
    for (i=0; i<SIDES_OF_ELEM(e); i++)
    {
      if (INNER_SIDE(e,i)) continue;

      /* allok a new boundary side */
      bs = geo->bs+geo->nbs;
      geo->nbs++;

      /* fill boundary side structure */
      bs->side = i;
      bs->nc = CORNERS_OF_SIDE(e,i);
      for (k=0; k<bs->nc; k++) bs->corners[k] = CORNER_OF_SIDE(e,i,k);

      /* select a boundary quadrature rule */
      if (DIM==2) {
        q = GetQuadrature(1,2,2);
        refFaceArea = 1.0;
      }
      if (DIM==3)
        switch (bs->nc)
        {
        case 3 :
          q = GetQuadrature(2,3,3);
          refFaceArea = 0.5;
          break;
        case 4 :
          q = GetQuadrature(2,4,4);
          refFaceArea = 1.0;
          break;
        default :
          return(1);
        }
      if (q==NULL) return(1);
      bs->nbgp = Q_NIP(q);

      /* collect global coordinates of corners of side */
      for (k=0; k<bs->nc; k++) {
        CP(sco_global[k],geo->co_global[bs->corners[k]]);
      }

      /* fill boundary gauss points */
      for (k=0; k<bs->nbgp; k++)
      {
        bgp = bs->bgp+k;
        bgp->weight = refFaceArea*Q_WEIGHT(q,k);
        for (j=0; j<DIM-1; j++) bgp->local[j] = Q_LOCAL(q,k)[j];

        /* interpolate param */
                                #ifdef __TWODIM__
        nvals[0] = 0.0;
        nvals[1] = 1.0;
        InterpolateFEFunction(DIM-1,bs->nc,bgp->local,nvals,&(bgp->param[0]));
                                #endif
                                #ifdef __THREEDIM__
        switch (bs->nc)
        {
        case 3 :
          nvals[0] = 0.0;
          nvals[1] = 1.0;
          nvals[2] = 0.0;
          InterpolateFEFunction(DIM-1,bs->nc,bgp->local,nvals,&(bgp->param[0]));
          nvals[0] = 0.0;
          nvals[1] = 0.0;
          nvals[2] = 1.0;
          InterpolateFEFunction(DIM-1,bs->nc,bgp->local,nvals,&(bgp->param[1]));                                  break;
        case 4 :
          nvals[0] = 0.0;
          nvals[1] = 1.0;
          nvals[2] = 1.0;
          nvals[3] = 0.0;
          InterpolateFEFunction(DIM-1,bs->nc,bgp->local,nvals,&(bgp->param[0]));
          nvals[0] = 0.0;
          nvals[1] = 0.0;
          nvals[2] = 1.0;
          nvals[3] = 1.0;
          InterpolateFEFunction(DIM-1,bs->nc,bgp->local,nvals,&(bgp->param[1]));                                  break;
          break;
        }
                                #endif

        /* calculate surface element */
        SurfaceElement(DIM,bs->nc,sco_global,bgp->local,&bgp->surfel);

        /* value of basis function on side at bgp */
        for (j=0; j<bs->nc; j++) nvals[j] = 0.0;
        for (j=0; j<bs->nc; j++) {
          nvals[j] = 1.0;
          InterpolateFEFunction(DIM-1,bs->nc,bgp->local,nvals,&(bgp->N[j]));
          nvals[j] = 0.0;
        }
      }
    }

  return(0);
}
