// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      fvgeom.h                                                          */
/*                                                                          */
/* Purpose:   geometry evaluation for Finite-Volume discretization			*/
/*			  dimension independent, general elmement						*/
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
/* History:   06.05.96 begin, ug version 3.2								*/
/*            06.08.96 modification and extension (Henrik Rentz-Reichert)   */
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

#ifndef __FVGEOM__
#define __FVGEOM__

#ifndef __GM__
#include "gm.h"
#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#ifndef MAXNC
#define MAXNC           MAX_CORNERS_OF_ELEM             /* just to make it more readable*/
#endif

#ifndef MAXE
#define MAXE            MAX_EDGES_OF_ELEM               /* just to make it more readable*/
#endif

#ifndef MAXS
#define MAXS            MAX_SIDES_OF_ELEM               /* just to make it more readable*/
#endif

#define MAXF            MAX_EDGES_OF_ELEM               /* max # scvf					*/
#define MAXBF           (MAX_SIDES_OF_ELEM*MAX_CORNERS_OF_SIDE)
/* max # boundary faces			*/

/* defines to specify data to be filled by EvaluateShapesAndDerivatives */
#define FILL_CORNER_DATA                        (1<<0)
#define FILL_SHAPES                                     (1<<1)
#define FILL_DERIV                                      (1<<2)
#define FILL_GRAD                                       (1<<3)
#define FILL_J                                      (1<<4)
#define FILL_SHAPES_AND_GRAD            (FILL_SHAPES | FILL_GRAD)
#define FILL_SHAPES_GRAD_AND_DERIV      (FILL_SHAPES | FILL_GRAD | FILL_DERIV)
#define FILL_ALL                                        (FILL_CORNER_DATA | FILL_SHAPES | FILL_GRAD | FILL_DERIV | FILL_J)

/* macros for access to exported data structures */
#define FVG_ELEM(p)                                     ((p)->e)
#define FVG_TAG(p)                                      ((p)->tag)
#define FVG_NSCV(p)                                     ((p)->n_scv)
#define FVG_NSCVF(p)                            ((p)->n_scvf)
#define FVG_NSCVBF(p)                           ((p)->n_bf)
#define FVG_GCO(p,i)                            ((p)->co_global[i])
#define FVG_GCOPTR(p)                           ((p)->co_global)
#define FVG_LCO(p,i)                            ((p)->co_local[i])
#define FVG_GEM(p,i)                            ((p)->em_global[i])
#define FVG_LEM(p,i)                            ((p)->em_local[i])
#define FVG_GSM(p,i)                            ((p)->sm_global[i])
#define FVG_LSM(p,i)                            ((p)->sm_local[i])
#define FVG_GCM(p)                                      ((p)->s_global)
#define FVG_LCM(p)                                      ((p)->s_local)
#define FVG_SCV(p,i)                            ((p)->scv+(i))
#define FVG_SCVF(p,i)                           ((p)->scvf+(i))
#define FVG_SCVBF(p,i)                          ((p)->bf+(i))
#define FVG_COSDV(p,i)                          ((p)->co_sdv+(i))
#define FVG_IPSDV(p,i)                          (&(((p)->scvf[i]).sdv))
#define FVG_BIPSDV(p,i)                         (&(((p)->bf[i]).sdv))

#define SCV_CO(p)                                       ((p)->co)
#define SCV_GCO(p)                                      ((p)->center)
#define SCV_VOL(p)                                      ((p)->volume)
#define SCV_NDPROP(p)                           ((p)->node_property)

#define SCVF_FROM(p)                            ((p)->i)
#define SCVF_TO(p)                                      ((p)->j)
#define SCVF_LIP(p)                                     ((p)->ip_local)
#define SCVF_GIP(p)                                     ((p)->ip_global)
#define SCVF_NORMAL(p)                          ((p)->normal)
#define SCVF_SDV(p)                                 ((p)->sdv)

#define SCVBF_FROM(p)                           ((p)->co)
#define SCVBF_SIDE(p)                           ((p)->side)
#define SCVBF_LIP(p)                            ((p)->ip_local)
#define SCVBF_PARAMPTR(p)                       ((p)->param)
#define SCVBF_PARAM(p,i)                        ((p)->param[i])
#define SCVBF_NORMAL(p)                         ((p)->normal)
#define SCVBF_AREA(p)                           ((p)->area)
#define SCVBF_SDV(p)                            ((p)->sdv)

#define SDV_SHAPEPTR(p)                         ((p)->shape)
#define SDV_SHAPE(p,i)                          ((p)->shape[i])
#define SDV_GRADPTR(p,i)                        ((p)->grad[i])
#define SDV_GRAD(p,i,j)                         ((p)->grad[i][j])
#define SDV_JINV(p)                                     ((p)->Jinv)
#define SDV_J(p)                                        ((p)->J)
#define SDV_DETJ(p)                                     ((p)->detJ)

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

/* shape functions and derivatives etc. */
typedef struct {
  DOUBLE shape[MAXNC];                          /* values of shape functions at ip		*/
  DOUBLE_VECTOR grad[MAXNC];                    /* derivatives of shape functions at ip	*/
  DOUBLE J[DIM][DIM];                               /* jacobian at ip			            */
  DOUBLE Jinv[DIM][DIM];                        /* inverse of jacobian at ip			*/
  DOUBLE detJ;                                          /* det of jacobian at ip				*/

} SD_VALUES;

/* geometry related data */
typedef struct {
  INT co;                                                       /* # of corner							*/
  DOUBLE_VECTOR center;                         /* node position						*/
  DOUBLE volume;                                        /* volume (area) of scv					*/
  INT node_property;                                    /* subdomain info						*/
} SubControlVolume;                                     /* FV intersected with element			*/

typedef struct {
  INT i,j;                                                      /* scvf seperates corner i and j of elem*/
  DOUBLE_VECTOR ip_local;                       /* integration point in local coords	*/
  DOUBLE_VECTOR ip_global;                          /* integration point in global coords	*/
  DOUBLE_VECTOR normal;                         /* normal on face at ip pointing to CV j*/
  SD_VALUES sdv;                                        /* shape fcts, deriv. etc. at scv-faces	*/
} SubControlVolumeFace;

typedef struct {
  INT co;                                                       /* corresponding corner					*/
  INT side;                                                     /* boundary side of element				*/
  DOUBLE_VECTOR ip_local;                       /* integration point in local coords	*/
  DOUBLE param[DIM-1];                          /* local side coordinates                       */
  DOUBLE_VECTOR normal;                         /* normal on face at ip pointing to CV j*/
  DOUBLE area;                                          /* area of boundary face				*/
  SD_VALUES sdv;                                        /* shape fcts, deriv. etc. at b-faces	*/
} BoundaryFace;

typedef struct {
  const ELEMENT *e;                                     /* data for this element				*/
  INT tag;                                                      /* element type							*/
  INT n_scv;                                                    /* # sub control volumes (==corners)	*/
  INT n_scvf;                                                   /* # sub control volume faces (==ip's)  */
  INT n_bf;                                                     /* # boundary faces						*/

  DOUBLE_VECTOR co_global[MAXNC];       /* points in global space, corners      */
  DOUBLE_VECTOR co_local[MAXNC];        /* points in local space, corners       */
  DOUBLE_VECTOR em_global[MAXE];        /* points in global space, edge midpoint*/
  DOUBLE_VECTOR em_local[MAXE];         /* points in local space, edge midpoints*/
  DOUBLE_VECTOR sm_global[MAXS];        /* points in global space, side midpoint*/
  DOUBLE_VECTOR sm_local[MAXS];         /* points in local space, side midpoints*/
  DOUBLE_VECTOR s_global;                       /* points in global space, center       */
  DOUBLE_VECTOR s_local;                        /* points in local space, center        */
  SD_VALUES co_sdv[MAXNC];                      /* shape fcts, deriv. etc. at corners	*/
  SubControlVolume scv[MAXNC];          /* sub control volumes					*/
  SubControlVolumeFace scvf[MAXF];      /* sub control volume faces				*/
  BoundaryFace bf[MAXBF];                       /* boundary faces						*/

} FVElementGeometry;                            /* geometry data for a general element	*/


/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

/* finite volume geometry */
INT EvaluateFVGeometry                          (const ELEMENT *e, FVElementGeometry *geo);

/* upwinding procedures */
INT GetFullUpwindShapes                         (const FVElementGeometry *geo, const DOUBLE_VECTOR IPVel[MAXF], DOUBLE Shape[MAXF][MAXNC]);
INT GetSkewedUpwindShapes                       (const FVElementGeometry *geo, const DOUBLE_VECTOR IPVel[MAXF], DOUBLE Shape[MAXF][MAXNC]);
INT GetMJRawRegularUpwindShapes         (const FVElementGeometry *geo, const DOUBLE_VECTOR IPVel[MAXF],
                                         DOUBLE NodalShape[MAXF][MAXNC], DOUBLE IPShape[MAXF][MAXF]);
INT GetMJRawPositiveUpwindShapes        (const FVElementGeometry *geo, const DOUBLE_VECTOR IPVel[MAXF],
                                         DOUBLE NodalShape[MAXF][MAXNC], DOUBLE IPShape[MAXF][MAXF]);

/* aligned finite volumes */
INT AFVGeometry                                         (const ELEMENT *theElement, FVElementGeometry *geo, DOUBLE_VECTOR Convection);

/* shape functions and their derivatives */
INT EvaluateShapesAndDerivatives        (FVElementGeometry *geo, INT flags);

/* intersect polygon with line */
INT Intersect2d (INT nco, const DOUBLE_VECTOR *x, const DOUBLE_VECTOR vel, const DOUBLE_VECTOR pt,
                 INT *Side, DOUBLE lambda[DIM_OF_BND]);

/* init */
INT InitFiniteVolumeGeom                        (void);

#endif
