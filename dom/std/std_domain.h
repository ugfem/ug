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


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __STD_DOMAIN__
#define __STD_DOMAIN__

#include "compiler.h"

#include "ugenv.h"

#include "domain.h"


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

typedef DOUBLE COORD_BND_VECTOR[DIM_OF_BND];

#define BVP_STANDARD                     0
#define BVP_MARC                         1

#define POINT_PATCH_TYPE                 0
#define LINE_PATCH_TYPE                  1
#define LINEAR_PATCH_TYPE                2
#define PARAMETRIC_PATCH_TYPE            3
#define MARC_0_PATCH_TYPE                4
#define MARC_1_PATCH_TYPE                5
#define MARC_2_PATCH_TYPE                6

/* macros for DOMAIN_PART_INFO */
#define DPI_SD2P_PTR(p)                                 ((p)->sd2part)
#define DPI_SD2P(p,sd)                                  ((p)->sd2part[sd])
#define DPI_SG2P_PTR(p)                                 ((p)->sg2part)
#define DPI_SG2P(p,sg)                                  ((p)->sg2part[sg])
#define DPI_LN2P_PTR(p)                                 ((p)->ln2part)
#define DPI_LN2P(p,c0,c1)                               ((p)->ln2part[c0][c1])
#define DPI_PT2P_PTR(p)                                 ((p)->pt2part)
#define DPI_PT2P(p,pt)                                  ((p)->pt2part[pt])

/* macros for DOMAIN */
#define DOMAIN_MIDPOINT(p)                              ((p)->MidPoint)
#define DOMAIN_RADIUS(p)                                ((p)->radius)
#define DOMAIN_CONVEX(p)                                ((p)->domConvex)
#define DOMAIN_NSEGMENT(p)                              ((p)->numOfSegments)
#define DOMAIN_NCORNER(p)                               ((p)->numOfCorners)
#define DOMAIN_NPARTS(p)                                ((p)->nParts)
#define DOMAIN_PARTINFO(p)                              ((p)->dpi)

/* macros for STD_BVP */
#define MAX_CORNERS_OF_LINEAR_PATCH      DIM
#undef  CORNERS_OF_BND_SEG
#define CORNERS_OF_BND_SEG               2*DIM_OF_BND

#define STD_BVP_NAME(p)                                 ENVITEM_NAME(p)

#define STD_BVP_DOMAIN(p)                               ((p)->Domain)
#define STD_BVP_PROBLEM(p)                              ((p)->Problem)

#define STD_BVP_MIDPOINT(p)                             ((p)->MidPoint)
#define STD_BVP_RADIUS(p)                               ((p)->radius)
#define STD_BVP_CONVEX(p)                               ((p)->domConvex)

#define STD_BVP_NCORNER(p)                              ((p)->ncorners)
#define STD_BVP_NSIDES(p)                               ((p)->nsides)
#define STD_BVP_SIDEOFFSET(p)                   ((p)->sideoffset)
#define STD_BVP_PATCHES(p)                              ((p)->patches)
#define STD_BVP_PATCH(p,i)                              ((p)->patches[i])
#define STD_BVP_CORNERPATCH(p,i)                ((p)->patches[i])
#define STD_BVP_LINEPATCH(p,i)                  ((p)->patches[i+STD_BVP_NCORNER(p)])
#define STD_BVP_SIDEPATCH(p,i)                  ((p)->patches[i+STD_BVP_SIDEOFFSET(p)])

#define STD_BVP_NSUBDOM(p)                              ((p)->numOfSubdomains)
#define STD_BVP_NDOMPART(p)                             ((p)->nDomainParts)
#define STD_BVP_S2P_PTR(p)                              ((p)->s2p)
#define STD_BVP_S2P(p,s)                                ((p)->s2p[s])
#define STD_BVP_NPATCH(p)                               ((p)->nPatch)
#define STD_BVP_NCOEFFPROC(p)                   ((p)->numOfCoeffFct)
#define STD_BVP_NUSERPROC(p)                    ((p)->numOfUserFct)
#define STD_BVP_CONFIGPROC(p)                   ((p)->ConfigProc)
#define STD_BVP_COEFFPROC(p,i)                  ((CoeffProcPtr)((p)->CU_ProcPtr[i]))
#define STD_BVP_USERPROC(p,i)                   ((UserProcPtr)((p)->CU_ProcPtr[i+STD_BVP_NCOEFFPROC(p)]))

#define GetSTD_BVP(p)                           ((STD_BVP *)(p))

/****************************************************************************/
/*																			*/
/* macros for boundary segments												*/
/*																			*/
/****************************************************************************/

#define SEG_LEFT(p)                                             ((p)->left)
#define SEG_RIGHT(p)                                    ((p)->right)
#define SEG_ID(p)                                               ((p)->id)
#define SEG_TYPE(p)                                             ((p)->segType)
#define SEG_POINTS(p)                                   ((p)->points)
#define SEG_POINT(p,i)                                  ((p)->points[i])
#define SEG_RESOLUTION(p)                               ((p)->resolution)
#define SEG_ALPHA(p)                                    ((p)->alpha)
#define SEG_FROM(p,i)                                   ((p)->alpha[i])
#define SEG_BETA(p)                                             ((p)->beta)
#define SEG_TO(p,i)                                             ((p)->beta[i])
#define SEG_FUNC(p)                                             ((p)->BndSegFunc)
#define SEG_DATA(p)                                             ((p)->data)

/****************************************************************************/
/*																			*/
/* macros for patches														*/
/*																			*/
/****************************************************************************/

#define PATCH_FIXED                             0
#define PATCH_BND_OF_FREE               1
#define PATCH_FREE                              2

#define PATCH_TYPE(p)           (p)->ge.type
#define PATCH_STATE(p)          (p)->ge.state
#define PATCH_IS_FREE(p)                ((p)->ge.state==PATCH_FREE)
#define PATCH_IS_FIXED(p)               ((p)->ge.state==PATCH_FIXED)
#define PATCH_ID(p)             (p)->ge.id
#define POINT_PATCH_N(p)        (p)->po.npatches
#define POINT_PATCH_PID(p,i)    (p)->po.pop[i].patch_id
#define POINT_PATCH_CID(p,i)    (p)->po.pop[i].corner_id
#define LINE_PATCH_C0(p)                ((p)->li.c0)
#define LINE_PATCH_C1(p)                ((p)->li.c1)
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
#define LINEAR_PATCH_LEFT(p)    (p)->lp.left
#define LINEAR_PATCH_RIGHT(p)   (p)->lp.right
#define LINEAR_PATCH_N(p)       (p)->lp.corners
#define LINEAR_PATCH_POINTS(p,i) (p)->lp.points[i]
#define LINEAR_PATCH_POS(p,i)    (p)->lp.pos[i]

#define BND_PATCH_ID(p)         (((BND_PS *)p)->patch_id)
#define BND_DATA(p)             (((BND_PS *)p)->data)
#define BND_N(p)                (((BND_PS *)p)->n)
#define BND_LOCAL(p,i)          (((BND_PS *)p)->local[i])
#define BND_SIZE(p)             ((((BND_PS *)p)->n-1)*sizeof(COORD_BND_VECTOR)+sizeof(BND_PS))
#define M_BNDS_NSIZE(n)         (((n)-1)*sizeof(M_BNDP)+sizeof(M_BNDS))
#define M_BNDS_SIZE(p)          M_BNDS_NSIZE(((M_BNDS *)(p))->n)


#define IF_MARC(p) \
  if (PATCH_TYPE(currBVP->patches[BND_PATCH_ID(p)]) >= MARC_0_PATCH_TYPE)

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

typedef INT (*BndSegFuncPtr)(void *,DOUBLE *,DOUBLE *);


/*----------- definition of structs ----------------------------------------*/

typedef struct {

  const INT *sd2part;                                           /* table subdomain to part			*/
  const INT *sg2part;                                           /* table segment   to part			*/
#       ifdef __THREEDIM__
  const INT **ln2part;                                  /* table line	   to part			*/
#       endif
  const INT *pt2part;                                           /* table point	   to part			*/

} DOMAIN_PART_INFO;

struct domain {

  /* fields for environment directory */
  ENVDIR d;

  /* domain variables */
  DOUBLE MidPoint[DIM];                                 /* point in the middle of domain	*/
  DOUBLE radius;                                                /* defines sphere around MidPoint	*/
  /* containing the domain			*/
  INT numOfSegments;                                            /* number of boundary segments		*/
  INT numOfCorners;                                             /* number of corner points			*/
  INT domConvex;                                                /* is the domain convex?			*/

  /* description of domain parts */
  INT nParts;                                                           /* number of parts in the domain	*/
  const DOMAIN_PART_INFO *dpi;                  /* domain part info					*/
};

struct boundary_segment {

  /* fields for environment directory */
  ENVVAR v;

  /* fields for boundary segment */
  INT left,right;                                         /* number of left and right subdomain */
  INT id;                                                         /* unique id of that segment			*/
  INT segType;                                            /* segment type, see above			*/
  INT points[CORNERS_OF_BND_SEG];         /* numbers of the vertices (ID)		*/
  INT resolution;                                         /* measure for the curvature			*/
  DOUBLE alpha[DIM_OF_BND],beta[DIM_OF_BND];              /* parameter interval used*/
  BndSegFuncPtr BndSegFunc;                       /* pointer to definition function     */
  void *data;                                             /* can be used by applic to find data */
};

struct linear_segment {

  /* fields for environment directory */
  ENVVAR v;

  /* fields for boundary segment */
  INT left,right;                                         /* number of left and right subdomain */
  INT id;                                                         /* unique id of that segment			*/
  INT n;                                  /* number of corners                  */
  INT points[MAX_CORNERS_OF_LINEAR_PATCH];       /* numbers of the vertices (ID)*/
  DOUBLE x[MAX_CORNERS_OF_LINEAR_PATCH][DIM_OF_BND];            /* coordinates  */
};

/****************************************************************************/
/*																			*/
/* problem data structure													*/
/*																			*/
/****************************************************************************/

/*----------- typedef for functions ----------------------------------------*/

typedef INT (*BndCondProcPtr)(void *, void *, DOUBLE *, DOUBLE *, INT *);

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
  INT type;
  struct domain *Domain;             /* domain pointer                      */
  struct problem *Problem;           /* problem pointer                     */

  char bnd_file[NAMESIZE];               /* file name for boundary infos        */
  char mesh_file[NAMESIZE];              /* file name for meshinfos             */

  /* domain part */
  DOUBLE MidPoint[DIM];              /* sphere in which the domain lies     */
  DOUBLE radius;
  INT domConvex;                     /* 1 if domain is convex, 0 if not     */
  INT numOfSubdomains;               /* nb. of subdomains,
                                                                                exterior not counted                */
  INT nDomainParts;                                      /* number of parts in the domain		*/
  INT *s2p;                                                      /* pointer to table subbdom --> part	*/

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
  INT state;                        /* fixed/bnd of free/free               */
  INT id;                           /* unique id used for load/store        */
};

struct point_on_patch {

  INT patch_id;
  INT corner_id;
};

struct point_patch {

  INT type;                         /* patch type                           */
  INT state;                        /* fixed/bnd of free/free               */
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
  INT state;                        /* fixed/bnd of free/free               */
  INT id;                           /* unique id used for load/store        */

  INT npatches;                     /* number of patches                    */
  INT c0;                                                       /* corner 0 of line						*/
  INT c1;                                                       /* corner 1 of line						*/
  struct line_on_patch lop[1];      /* reference to surface                 */
};
#endif

struct linear_patch {

  INT type;                         /* patch type                           */
  INT state;                        /* fixed/bnd of free/free               */
  INT id;                           /* unique id used for load/store        */

  INT left,right;                                       /* id of left and right subdomain       */
  INT corners;                      /* number of corners                    */
  INT points[MAX_CORNERS_OF_LINEAR_PATCH];    /* ids of points              */
  DOUBLE pos[MAX_CORNERS_OF_LINEAR_PATCH][DIM];   /* position               */
};

struct parameter_patch {

  INT type;                         /* patch type                           */
  INT state;                        /* fixed/bnd of free/free               */
  INT id;                           /* unique id used for load/store        */

  INT left,right;                                       /* id of left and right subdomain       */
  INT points[CORNERS_OF_BND_SEG];   /* ids of points                        */
  DOUBLE range[2][DIM_OF_BND];       /* parameter range                      */

  BndSegFuncPtr BndSegFunc;                     /* pointer to definition function           */
  void *bs_data;                                        /* can be used by applic to find data   */

  /* fields for boundary condition */
  BndCondProcPtr BndCond;                   /* function defining boundary condition */
  void *bc_data;                                    /* additional data for bnd cond             */
};

struct marc_0_patch {

  INT type;                         /* patch type                           */
  INT state;                        /* fixed/bnd of free/free               */
  INT id;                           /* unique id used for load/store        */

  DOUBLE pos[3];                        /* position                             */
};

struct marc_1_patch {

  INT type;                         /* patch type                           */
  INT state;                        /* fixed/bnd of free/free               */
  INT id;                           /* unique id used for load/store        */

  INT p[2];                             /* line between two points              */
};

struct marc_2_patch {

  INT type;                         /* patch type                           */
  INT state;                        /* fixed/bnd of free/free               */
  INT id;                           /* unique id used for load/store        */

  INT c;                                /* bnd cond                             */
  INT p[3];                             /* triangle of three points             */
};

union patch {
  struct generic_patch ge;
  struct point_patch po;
  struct linear_patch lp;
  struct parameter_patch pa;
    #ifdef __THREEDIM__
  struct line_patch li;
        #endif
  struct marc_0_patch m0;
  struct marc_1_patch m1;
  struct marc_2_patch m2;
} ;

struct bnd_ps {

  INT patch_id;                     /* associated patch                     */
  void *data;                                           /* e.g. global coordiantes, pointers...	*/
  INT n;                            /* number of arguments                  */
  COORD_BND_VECTOR local[1];        /* parameter range                      */
};

struct marc_bndp {

  INT patch_id;                     /* associated patch                     */
  DOUBLE pos[3];                                /* position                             */
};

struct marc_bnds {

  INT patch_id;                     /* associated patch                     */
  INT n;                                /* number of corners                    */
  struct marc_bndp *p[1];               /* corners                              */
};

/*----------- typedef for structs ------------------------------------------*/

/* typedefs */
#undef DOMAIN
typedef struct domain DOMAIN;
typedef struct linear_segment LINEAR_SEGMENT;
typedef struct boundary_segment BOUNDARY_SEGMENT;
typedef struct problem PROBLEM;
typedef struct bndcond BOUNDARY_CONDITION;

typedef struct std_BoundaryValueProblem STD_BVP;
typedef union patch PATCH;
typedef struct point_patch POINT_PATCH;
typedef struct line_patch LINE_PATCH;
typedef struct linear_patch LINEAR_PATCH;
typedef struct parameter_patch PARAMETER_PATCH;
typedef struct marc_0_patch M0_PATCH;
typedef struct marc_1_patch M1_PATCH;
typedef struct marc_2_patch M2_PATCH;
typedef struct bnd_ps BND_PS;
typedef struct marc_bndp M_BNDP;
typedef struct marc_bnds M_BNDS;

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

#       ifdef __THREEDIM__
INT RepairMesh (HEAP *Heap, INT MarkKey, MESH *mesh);
#       endif


/* domain definition */
DOMAIN                     *CreateDomainWithParts       (char *name, DOUBLE *MidPoint,
                                                         DOUBLE radius, INT segments,
                                                         INT corners, INT Convex,
                                                         INT nParts, const DOMAIN_PART_INFO *dpi);
DOMAIN                     *CreateDomain                        (char *name, DOUBLE *MidPoint,
                                                                 DOUBLE radius, INT segments,
                                                                 INT corners, INT Convex);
DOMAIN *GetDomain                           (char *name);
BOUNDARY_SEGMENT   *CreateBoundarySegment       (char *name, INT left, INT right,
                                                 INT id, INT type, INT res,
                                                 INT *point,
                                                 DOUBLE *alpha, DOUBLE *beta,
                                                 BndSegFuncPtr BndSegFunc,
                                                 void *data);
BOUNDARY_SEGMENT   *CreateBoundarySegment2D     (char *name, int left, int right,
                                                 int id, int from, int to, int res,
                                                 DOUBLE alpha, DOUBLE beta,
                                                 BndSegFuncPtr BndSegFunc,
                                                 void *data);
LINEAR_SEGMENT *CreateLinearSegment (char *name,
                                     INT left, INT right,INT id,
                                     INT n, INT *point,
                                     DOUBLE x[MAX_CORNERS_OF_LINEAR_PATCH][DIM]);

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
BVP   *Create_MarcBVP             (char *BVPname, BndCondProcPtr theBndCond,
                                   int numOfCoeffFct, CoeffProcPtr coeffs[],
                                   int numOfUserFct, UserProcPtr userfct[]);

/* scanning of coordinates */
INT   ReadAndPrintArgvPosition    (char *name, INT argc, char **argv, DOUBLE *pos);

#endif
