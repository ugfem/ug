// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  std_domain.h													*/
/*																			*/
/* Purpose:   standard domain declaration	                                                                */
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

typedef COORD COORD_BND_VECTOR[DIM_OF_BND];

#define POINT_PATCH_TYPE                 0
#define LINE_PATCH_TYPE                  1
#define LINEAR_PATCH_TYPE                2
#define PARAMETRIC_PATCH_TYPE            3

#define MAX_CORNERS_OF_LINEAR_PATCH      4
#undef  CORNERS_OF_BND_SEG
#define CORNERS_OF_BND_SEG               2*DIM_OF_BND

/****************************************************************************/
/*																			*/
/* macros for patches														*/
/*																			*/
/****************************************************************************/

#define PATCH_TYPE(p)           (p)->ge.type
#define PATCH_ID(p)             (p)->ge.id
#define POINT_PATCH_N(p)        (p)->po.npatches
#define POINT_PATCH_PID(p,i)    (p)->po.pop[i].patch_id
#define POINT_PATCH_CID(p,i)    (p)->po.pop[i].corner_id
#define LINE_PATCH_N(p)         (p)->li.npatches
#define LINE_PATCH_PID(p,i)     (p)->li.lop[i].patch_id
#define LINE_PATCH_CID0(p,i)    (p)->li.lop[i].corner_id[0]
#define LINE_PATCH_CID1(p,i)    (p)->li.lop[i].corner_id[1]
#define PARAM_PATCH_LEFT(p)     (p)->pa.left
#define PARAM_PATCH_RIGHT(p)    (p)->pa.right
#define PARAM_PATCH_POINTS(p,i) (p)->pa.points[i]
#define PARAM_PATCH_RANGE(p)    (p)->pa.range
#define PARAM_PATCH_BS(p)       (p)->pa.BndSegFunc
#define PARAM_PATCH_BSD(p)      (p)->pa.bs_data
#define PARAM_PATCH_BC(p)       (p)->pa.BndCond
#define PARAM_PATCH_BCD(p)      (p)->pa.bc_data

#define BND_PATCH_ID(p)         ((BND_PS *)p)->patch_id
#define BND_N(p)                ((BND_PS *)p)->n
#define BND_LOCAL(p,i)          ((BND_PS *)p)->local[i]
#define BND_SIZE(p)             ((((BND_PS *)p)->n-1)*sizeof(COORD_BND_VECTOR)+sizeof(BND_PS))

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

typedef INT (*BndCondProcPtr)(void *, void *, COORD *, DOUBLE *, INT *);

/*----------- definition of structs ----------------------------------------*/

struct problem {

  /* fields for environment directory */
  ENVDIR d;

  /* fields for problem */
  INT problemID;                                /* used to identify problem type			*/
  ConfigProcPtr ConfigProblem;      /* procedure to reinitialize problem		*/
  INT numOfCoeffFct;                            /* number of coefficient functions			*/
  INT numOfUserFct;                             /* number of User functions					*/
  void *CU_ProcPtr[1];                  /* coefficient functions					*/
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
/* BoundaryValueProblem data structure										*/
/*																			*/
/****************************************************************************/

struct std_BoundaryValueProblem
{
  /* fields for environment directory */
  ENVDIR d;

  /* init */
  struct domain *Domain;             /* domain pointer                      */
  struct problem *Problem;           /* problem pointer                     */

  /* domain part */
  COORD MidPoint[DIM];               /* sphere in which the domain lies     */
  COORD radius;
  INT domConvex;                     /* 1 if domain is convex, 0 if not     */
  INT numOfSubdomains;               /* nb. of subdomains,
                                                                                exterior not counted                */
  /* boundary decription */
  INT ncorners;
  INT nsides;
  INT sideoffset;
  union patch **patches;                     /* list of patches	                            */

  /* problem part */
  ConfigProcPtr ConfigProc;          /* configuration function              */
  INT numOfCoeffFct;                 /* nb. of coefficient functions        */
  INT numOfUserFct;                  /* nb. of user functions               */

  BndCondProcPtr GeneralBndCond;     /* general bnd. cond. (if exists)      */
  void *CU_ProcPtr[1];                       /* coefficient functions				*/
};

/****************************************************************************/
/*																			*/
/* Patch data structure														*/
/*																			*/
/****************************************************************************/

struct generic_patch {

  INT type;                         /* patch type                           */
  INT id;                           /* unique id used for load/store        */
};

struct point_on_patch {

  INT patch_id;
  INT corner_id;
};

struct point_patch {

  INT type;                         /* patch type                           */
  INT id;                           /* unique id used for load/store        */

  INT npatches;                     /* number of patches                    */
  struct point_on_patch pop[1];     /* reference to surface                 */
};

#ifdef __THREEDIM__
struct line_on_patch {

  INT patch_id;
  INT corner_id[2];
};

struct line_patch {

  INT type;                         /* patch type                           */
  INT id;                           /* unique id used for load/store        */

  INT npatches;                     /* number of patches                    */
  struct line_on_patch lop[1];      /* reference to surface                 */
};
#endif

struct linear_patch {

  INT type;                         /* patch type                           */
  INT id;                           /* unique id used for load/store        */

  INT corners;                      /* number of corners                    */
  COORD pos[MAX_CORNERS_OF_LINEAR_PATCH][DIM];   /* position              */

  /* fields for boundary condition */
  BndCondProcPtr BndCond;                   /* function defining boundary condition */
  void *bc_data;                                    /* additional data for bnd cond             */
};

struct parameter_patch {

  INT type;                         /* patch type                           */
  INT id;                           /* unique id used for load/store        */

  INT left,right;                                       /* id of left and right subdomain       */
  INT points[CORNERS_OF_BND_SEG];   /* ids of points                        */
  COORD range[2][DIM_OF_BND];       /* parameter range                      */

  BndSegFuncPtr BndSegFunc;                     /* pointer to definition function           */
  void *bs_data;                                        /* can be used by applic to find data   */

  /* fields for boundary condition */
  BndCondProcPtr BndCond;                   /* function defining boundary condition */
  void *bc_data;                                    /* additional data for bnd cond             */
};

struct bnd_ps {

  INT patch_id;                     /* associated patch                     */
  INT n;                            /* number of arguments                  */
  COORD_BND_VECTOR local[1];        /* parameter range                      */
};

union patch {
  struct generic_patch ge;
  struct point_patch po;
  struct parameter_patch pa;
    #ifdef __THREEDIM__
  struct line_patch li;
        #endif
} ;

/*----------- typedef for structs ------------------------------------------*/

/* typedefs */
#undef DOMAIN
typedef struct domain DOMAIN;
typedef struct boundary_segment BOUNDARY_SEGMENT;
typedef struct problem PROBLEM;
typedef struct bndcond BOUNDARY_CONDITION;

typedef struct std_BoundaryValueProblem STD_BVP;
typedef union patch PATCH;
typedef struct point_patch POINT_PATCH;
typedef struct line_patch LINE_PATCH;
typedef struct linear_patch LINEAR_PATCH;
typedef struct parameter_patch PARAMETER_PATCH;
typedef struct bnd_ps BND_PS;

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
DOMAIN                     *CreateDomain                        (char *name, COORD *MidPoint,
                                                                 COORD radius, INT segments,
                                                                 INT corners, INT Convex);
DOMAIN *GetDomain                           (char *name);
BOUNDARY_SEGMENT   *CreateBoundarySegment       (char *name, INT left, INT right,
                                                 INT id, INT type, INT res,
                                                 INT *point,
                                                 COORD *alpha, COORD *beta,
                                                 BndSegFuncPtr BndSegFunc,
                                                 void *data);
BOUNDARY_SEGMENT   *CreateBoundarySegment2D     (char *name, int left, int right,
                                                 int id, int from, int to, int res,
                                                 COORD alpha, COORD beta,
                                                 BndSegFuncPtr BndSegFunc,
                                                 void *data);

/* problem definition */
PROBLEM                    *CreateProblem                       (char *domain, char *name,
                                                                 int id, ConfigProcPtr config,
                                                                 int numOfCoefficients,
                                                                 CoeffProcPtr coeffs[],
                                                                 int numOfUserFct,
                                                                 UserProcPtr userfct[]);
BOUNDARY_CONDITION *CreateBoundaryCondition (char *name, INT id,
                                             BndCondProcPtr theBndCond,
                                             void *Data);

/* BVP definition */
INT   STD_BVP_Configure           (INT argc, char **argv);
BVP   *CreateBoundaryValueProblem (char *BVPname, BndCondProcPtr theBndCond,
                                   int numOfCoeffFct, CoeffProcPtr coeffs[],
                                   int numOfUserFct, UserProcPtr userfct[]);
BVP       *CreateBVP                              (char *BVP, char *Domain, char *Problem);

#endif
