// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  header.h														*/
/*																			*/
/* Purpose:   standard header file template                                                             */
/*																			*/
/* Author:	  Peter Bastian/Klaus Johannsen                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __STD_DOMAIN__
#define __STD_DOMAIN__

#ifndef __SWITCH__
#include "switch.h"
#endif

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __DOMAIN__
#include "domain.h"
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

/****************************************************************************/
/*																			*/
/* macros for BVP															*/
/*																			*/
/****************************************************************************/

#define PROBLEMID(p)   (p)->problemID
#define CONFIGFUNC(p)  (p)->ConfigProblem
#define NUMCOEFF(p)    (p)->numOfCoeffFct
#define NUMUSERFCT(p)  (p)->numOfUserFct
#define COEFFFUNC(p,i) (CoeffProcPtr)((p)->CU_ProcPtr[i])
#define USERFUNC(p,i)  (UserProcPtr)((p)->CU_ProcPtr[i+(p)->numOfCoeffFct])

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* domain definition data structures										*/
/*																			*/
/****************************************************************************/

/*----------- typedef for functions ----------------------------------------*/

typedef INT (*BndSegFuncPtr)(void *,COORD *,COORD *);


/*----------- definition of structs ----------------------------------------*/

struct domain {

  /* fields for environment directory */
  ENVDIR d;

  /* domain variables */
  COORD MidPoint[DIM];                                  /* point in the middle of domain	*/
  COORD radius;                                                 /* defines sphere around MidPoint	*/
  /* containing the domain			*/
  INT numOfSegments;                                            /* number of boundary segments		*/
  INT numOfCorners;                                             /* number of corner points			*/
  INT domConvex;                                                /* is the domain convex?			*/
} ;

struct boundary_segment {

  /* fields for environment directory */
  ENVVAR v;

  /* fields for boundary segment */
  INT left,right;                                         /* number of left and right subdomain */
  INT id;                                                         /* unique id of that segment			*/
  INT segType;                                            /* segment type, see above			*/
  INT points[CORNERS_OF_BND_SEG];         /* numbers of the vertices (ID)		*/
  INT resolution;                                         /* measure for the curvature			*/
  COORD alpha[DIM_OF_BND],beta[DIM_OF_BND];               /* parameter interval used*/
  BndSegFuncPtr BndSegFunc;                       /* pointer to definition function     */
  void *data;                                             /* can be used by applic to find data */
} ;

/****************************************************************************/
/*																			*/
/* problem data structure													*/
/*																			*/
/****************************************************************************/

/*----------- typedef for functions ----------------------------------------*/

#ifdef __version23__
typedef INT (*BndCondProcPtr)(void *, DOUBLE *, DOUBLE *, INT *);
#else
typedef INT (*BndCondProcPtr)(void *, void *, COORD *, DOUBLE *, INT *);
#endif

/*----------- definition of structs ----------------------------------------*/

struct problem {

  /* fields for environment directory */
  ENVDIR d;

  /* fields for problem */
  INT problemID;                                /* used to identify problem type			*/
  ConfigProcPtr ConfigProblem;      /* procedure to reinitialize problem		*/
  INT numOfCoeffFct;                            /* # of coefficient functions				*/
  INT numOfUserFct;                             /* # of User functions						*/
  void * CU_ProcPtr[1];                 /* coefficient functions					*/
};

struct bndcond {

  /* fields for environment variable */
  ENVVAR v;

  /* fields for boundary condition */
  INT id;                                               /* corresponds to boundary segment id !         */
  BndCondProcPtr BndCond;               /* function defining boundary condition         */
  void *data;                                   /* additional data for bnd cond                         */
};

/****************************************************************************/
/*																			*/
/* BondaryValueProblem data structure										*/
/*																			*/
/****************************************************************************/

struct std_BondaryValueProblem
{
  /* fields for environment directory */
  ENVDIR d;

  /* fields for domain */
  COORD MidPoint[DIM];                                  /* point in the middle of domain	*/
  COORD radius;                                                 /* defines sphere around MidPoint	*/
  /* containing the domain			*/
  INT numOfPatches;                                             /* number of patches				*/
  INT numOfCorners;                                             /* number of corner points			*/
  INT domConvex;                                                /* is the domain convex?			*/
  INT numOfSubdomains;                                  /* number of subdomains		        */

  /* fields for problem */
  INT problemID;                                /* used to identify problem type			*/
  ConfigProcPtr ConfigProblem;      /* procedure to reinitialize problem		*/
  INT numOfCoeffFct;                            /* # of coefficient functions				*/
  INT numOfUserFct;                             /* # of User functions						*/
  void * CU_ProcPtr[1];                 /* coefficient functions					*/
};

/****************************************************************************/
/*																			*/
/* Patch data structure														*/
/*																			*/
/****************************************************************************/

struct std_Patch {

  /* fields for environment directory */
  ENVVAR v;

  /* fields for boundary segment */
  INT left,right;                                         /* number of left and right subdomain */
  INT id;                                                         /* unique id of that patch			*/
  INT segType;                                            /* segment type, see above			*/
  INT nPoints;                                            /* numerber of corners of the segment */
  INT points[CORNERS_OF_BND_SEG];         /* numbers of the vertices (ID)		*/
  INT resolution;                                         /* measure for the curvature			*/
  COORD alpha[DIM_OF_BND],beta[DIM_OF_BND];               /* parameter interval used*/
  BndSegFuncPtr BndSegFunc;                       /* pointer to definition function     */
  void *bs_data;                                          /* can be used by applic to find data */

  /* fields for boundary condition */
  BndCondProcPtr BndCond;               /* function defining boundary condition         */
  void *bc_data;                                /* additional data for bnd cond                         */
} ;

/*----------- typedef for structs ------------------------------------------*/

/* typedefs */
typedef struct domain DOMAIN;
typedef struct boundary_segment BOUNDARY_SEGMENT;
typedef struct problem PROBLEM ;
typedef struct bndcond BOUNDARY_CONDITION;

typedef struct std_BondaryValueProblem STD_BVP;
typedef struct std_Patch STD_PATCH;

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* domain definition */
DOMAIN                   *CreateDomain                          (char *name, COORD *MidPoint, COORD radius, INT segments, INT corners, INT Convex);
DOMAIN                   *GetDomain                             (char *name);
BOUNDARY_SEGMENT *CreateBoundarySegment         (char *name, INT left, INT right,INT id,INT type,INT res,INT *point, COORD *alpha,COORD *beta, BndSegFuncPtr BndSegFunc, void *data);
BOUNDARY_SEGMENT *CreateBoundarySegment2D       (char *name, int left, int right, int id, int from, int to, int res, COORD alpha, COORD beta, BndSegFuncPtr BndSegFunc, void *data);
BOUNDARY_SEGMENT *GetFirstBoundarySegment       (DOMAIN *theDomain);
BOUNDARY_SEGMENT *GetNextBoundarySegment        (BOUNDARY_SEGMENT *theBSeg);
BOUNDARY_CONDITION *GetFirstBoundaryCondition (PROBLEM *theProblem);
BOUNDARY_CONDITION *GetNextBoundaryCondition (BOUNDARY_CONDITION *theBCond);

/* problem definition */
PROBLEM                  *CreateProblem                         (char *domain,char *name, int id, ConfigProcPtr config, int numOfCoefficients, CoeffProcPtr coeffs[], int numOfUserFct, UserProcPtr userfct[]);
PROBLEM                  *GetProblem                            (const char * domain, const char *name);
#ifdef __version23__
BOUNDARY_CONDITION *CreateBoundaryCondition (char *name, INT id, BndCondProcPtr theBndCond);
#endif
#ifdef __version3__
BOUNDARY_CONDITION *CreateBoundaryCondition (char *name, INT id, BndCondProcPtr theBndCond, void *Data);
#endif

/* BVP definition */
BVP                             *CreateBVP                                      (char *BVP, char *Domain, char *Problem);

#endif
