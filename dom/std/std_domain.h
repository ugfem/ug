// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \defgroup std The Standard Domain
 * \ingroup dom
 */
/*! \file std_domain.h
 * \ingroup std
 */

/** \addtogroup std
 *
 * @{
 */

/****************************************************************************/
/*                                                                          */
/* File:      std_domain.h                                                  */
/*                                                                          */
/* Purpose:   standard domain declaration                                   */
/*                                                                          */
/* Author:    Peter Bastian/Klaus Johannsen                                 */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70550 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   29.01.92 begin, ug version 2.0                                */
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

#ifndef __STD_DOMAIN__
#define __STD_DOMAIN__

#include "domain.h"

#include "dimension.h"
#include "namespace.h"

START_UGDIM_NAMESPACE

#undef  CORNERS_OF_BND_SEG
#define CORNERS_OF_BND_SEG               2*DIM_OF_BND


/** \todo Please doc me! */
typedef struct {

  /** \brief Table subdomain to part */
  const INT *sd2part;

  /** \brief Table segment to part */
  const INT *sg2part;

# ifdef __THREEDIM__
  /** \brief Table line to part */
  const INT **ln2part;
# endif

  /** \brief Table point to part */
  const INT *pt2part;

} DOMAIN_PART_INFO;

/*----------- typedef for functions ----------------------------------------*/
/** \brief ???
 *
 * \todo Please doc me!
 */
typedef INT (*BndSegFuncPtr)(void *,DOUBLE *,DOUBLE *);

/** \brief ???
 *
 * \todo Please doc me!
 */
typedef INT (*BndCondProcPtr)(void *, void *, DOUBLE *, DOUBLE *, INT *);


/* --- public functions --- */


/* domain definition */
void                     *CreateDomainWithParts       (const char *name, DOUBLE *MidPoint,
                                                       DOUBLE radius, INT segments,
                                                       INT corners, INT Convex,
                                                       INT nParts, const DOMAIN_PART_INFO *dpi);

void                     *CreateDomain                        (const char *name, DOUBLE *MidPoint,
                                                               DOUBLE radius, INT segments,
                                                               INT corners, INT Convex);

void   *CreateBoundarySegment       (char *name, INT left, INT right,
                                     INT id, enum BoundaryType type, INT res,
                                     INT *point,
                                     DOUBLE *alpha, DOUBLE *beta,
                                     BndSegFuncPtr BndSegFunc,
                                     void *data);

void   *CreateBoundarySegment2D     (char *name, int left, int right,
                                     int id, int from, int to, int res,
                                     DOUBLE alpha, DOUBLE beta,
                                     BndSegFuncPtr BndSegFunc,
                                     void *data);

/* problem definition */
void                    *CreateProblem                       (char *domain, char *name,
                                                              int id, ConfigProcPtr config,
                                                              int numOfCoefficients,
                                                              CoeffProcPtr coeffs[],
                                                              int numOfUserFct,
                                                              UserProcPtr userfct[]);

void *CreateBoundaryCondition (char *name, INT id,
                               BndCondProcPtr theBndCond,
                               void *Data);

BVP   *CreateBoundaryValueProblem (const char *BVPname, BndCondProcPtr theBndCond,
                                   int numOfCoeffFct, CoeffProcPtr coeffs[],
                                   int numOfUserFct, UserProcPtr userfct[]);
BVP       *CreateBVP                              (char *BVP, char *Domain, char *Problem);
const char *GetBVP_DomainName     (const BVP *aBVP);
const char *GetBVP_ProblemName    (const BVP *aBVP);
BVP   *Create_MarcBVP             (char *BVPname, BndCondProcPtr theBndCond,
                                   int numOfCoeffFct, CoeffProcPtr coeffs[],
                                   int numOfUserFct, UserProcPtr userfct[]);

BVP *CreateBVP_Problem (char *BVPName, char *DomainName, char *ProblemName);

END_NAMESPACE

/** @} */

#endif
